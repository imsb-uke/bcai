"""
Literature pipeline — synchronous orchestration functions.

No MCP, no async. Usable standalone or imported by the literature sub-agent server.

LLM model is read from the LLM_MODEL_LITERATURE env var (set in bcai/env):
    LLM_MODEL_LITERATURE = 'openrouter|openai/gpt-4o-mini'
Falls back to 'openai|gpt-4o-mini' if not set.

progress_cb: optional callable(str) passed in by the server for job progress updates.
"""

import os
import json
import time
import xml.etree.ElementTree as ET
from datetime import datetime

import httpx
import pandas as pd
from openai import OpenAI
from ollama import Client as OllamaClient

PUBMED_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
S2_SEARCH      = "https://api.semanticscholar.org/graph/v1/paper/search"


# ── LLM client (mirrors make_client in main/client.py, sync version) ─────────

def _make_client(llm_source: str):
    """Sync version of client.py's make_client. Returns (client, client_type)."""
    if llm_source == 'openai':
        return OpenAI(), 'openai'
    elif llm_source == 'ollama_local':
        return OpenAI(base_url="http://localhost:11434/v1", api_key="ollama"), 'openai'
    elif llm_source == 'ollama_cloud':
        return OllamaClient(
            host="https://ollama.com",
            headers={'Authorization': 'Bearer ' + os.getenv("OLLAMA_API_KEY", "")}
        ), 'ollama'
    elif llm_source == 'openrouter':
        return OpenAI(
            base_url="https://openrouter.ai/api/v1",
            api_key=os.getenv("OPENROUTER_API_KEY"),
        ), 'openai'
    else:
        raise ValueError(f"Unknown LLM source: '{llm_source}'. Use: openai, openrouter, ollama_local, ollama_cloud")


def _get_llm():
    """Parse LLM_MODEL_LITERATURE env var and return (client, client_type, model_name)."""
    model_str   = os.getenv("LLM_MODEL_LITERATURE", "openai|gpt-4o-mini")
    llm_source, model_name = model_str.split("|", 1)
    client, client_type    = _make_client(llm_source)
    return client, client_type, model_name


# ── Internal: PubMed fetcher ──────────────────────────────────────────────────

def _fetch_pubmed(query: str, max_results: int) -> list[dict]:
    resp = httpx.get(PUBMED_ESEARCH, params={
        "db": "pubmed", "term": query, "retmax": max_results, "retmode": "json",
    }, timeout=30)
    resp.raise_for_status()
    pmids = resp.json()["esearchresult"]["idlist"]
    if not pmids:
        return []

    resp = httpx.get(PUBMED_EFETCH, params={
        "db": "pubmed", "id": ",".join(pmids), "rettype": "xml", "retmode": "xml",
    }, timeout=60)
    resp.raise_for_status()

    papers = []
    for article in ET.fromstring(resp.text).findall(".//PubmedArticle"):
        pmid     = article.findtext(".//PMID", "")
        title    = article.findtext(".//ArticleTitle", "")
        abstract = " ".join(el.text or "" for el in article.findall(".//AbstractText"))
        year     = (article.findtext(".//PubDate/Year")
                    or article.findtext(".//PubDate/MedlineDate", "")[:4])
        authors  = [
            f"{a.findtext('ForeName', '')} {a.findtext('LastName', '')}".strip()
            for a in article.findall(".//Author")
        ]
        papers.append({
            "pmid": pmid, "title": title, "abstract": abstract,
            "year": year, "authors": authors[:5], "source": "pubmed",
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        })
    return papers


# ── Internal: Semantic Scholar fetcher ───────────────────────────────────────

def _fetch_semantic_scholar(query: str, max_results: int) -> list[dict]:
    headers = {"x-api-key": os.getenv("S2_API_KEY")} if os.getenv("S2_API_KEY") else {}  # claude
    # Retry up to 3 times with exponential back-off to handle 429 rate-limit responses
    for attempt in range(3):
        resp = httpx.get(S2_SEARCH, params={
            "query": query, "limit": min(max_results, 100),
            "fields": "title,abstract,year,authors,externalIds,citationCount",
        }, headers=headers, timeout=30)  # claude: added headers
        if resp.status_code == 429:
            time.sleep(5 * (3 ** attempt))   # claude: 5 s, 15 s, 45 s  (was 1, 2, 4)
            continue
        resp.raise_for_status()
        break
    else:
        resp.raise_for_status()   # all retries exhausted

    papers = []
    for p in resp.json().get("data", []):
        ext = p.get("externalIds") or {}
        papers.append({
            "pmid":           ext.get("PubMed", ""),
            "title":          p.get("title", ""),
            "abstract":       p.get("abstract") or "",
            "year":           str(p.get("year", "")),
            "authors":        [a["name"] for a in (p.get("authors") or [])[:5]],
            "citation_count": p.get("citationCount", 0),
            "source":         "semantic_scholar",
            "url":            f"https://www.semanticscholar.org/paper/{p.get('paperId', '')}",
        })
    return papers


# ── search_literature ─────────────────────────────────────────────────────────

def search_literature(
    query:       str,
    source:      str = "both",
    max_results: int = 20,
    file_dir:    str = "files",
    progress_cb       = None,
) -> dict:
    """
    Search PubMed and/or Semantic Scholar and save results as a JSON file.

    Args:
        query:       free-text search query (formulated by the mother LLM).
        source:      'pubmed' | 'semantic_scholar' | 'both' (default: 'both').
        max_results: max papers per source.
        file_dir:    directory to save output JSON.
        progress_cb: optional callable(str) for progress messages.

    Returns:
        {"n_papers": int, "papers_file": str, "next_steps": [...]}
    """
    def _p(msg):
        if progress_cb: progress_cb(msg)

    os.makedirs(file_dir, exist_ok=True)
    papers = []

    if source in ("pubmed", "both"):
        _p("searching PubMed...")
        papers += _fetch_pubmed(query, max_results)

    if source in ("semantic_scholar", "both"):
        _p("searching Semantic Scholar...")
        # claude ---
        try:
            papers += _fetch_semantic_scholar(query, max_results)
        except Exception as e:
            _p(f"Semantic Scholar unavailable ({e}); continuing with PubMed results only.")
        # ---

    # Deduplicate by PMID (or title for papers without one)
    seen, unique = set(), []
    for p in papers:
        key = p["pmid"] if p["pmid"] else p["title"]
        if key not in seen:
            seen.add(key)
            unique.append(p)

    safe_query = "".join(c if c.isalnum() else "_" for c in query)[:40]
    timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_file   = os.path.join(file_dir, f"literature_{safe_query}_{timestamp}.json")
    with open(out_file, "w") as f:
        json.dump(unique, f, indent=2)

    _p(f"done: {len(unique)} papers saved")
    return {
        "n_papers":    len(unique),
        "papers_file": out_file,
        "next_steps":  ["extract_information(papers_file, question)"],
    }


# ── Internal: single-paper LLM extraction ────────────────────────────────────

_EXTRACTION_PROMPT = """\
You are a biomedical research assistant. Read the paper abstract below and \
answer the following question as precisely as possible.

Question: {question}

Title: {title}
Abstract: {abstract}

Return a JSON object with exactly these fields:
  "answer"     – your answer to the question (null if the abstract does not contain relevant information)
  "evidence"   – the exact sentence(s) from the abstract that support your answer (null if none)
  "confidence" – "high", "medium", or "low"

Return ONLY the JSON object.\
"""


def _call_llm(client, client_type: str, model_name: str, prompt: str) -> dict:
    """Call the LLM and return parsed JSON. Handles both openai and ollama client types."""
    if client_type == 'openai':
        resp = client.chat.completions.create(
            model           = model_name,
            messages        = [{"role": "user", "content": prompt}],
            response_format = {"type": "json_object"},
            temperature     = 0,
        )
        return json.loads(resp.choices[0].message.content)
    else:  # ollama
        resp = client.chat(
            model    = model_name,
            messages = [{"role": "user", "content": prompt}],
            stream   = False,
        )
        raw = resp.message.content
        # Ollama doesn't enforce JSON format — try to parse, fall back to raw text
        try:
            return json.loads(raw)
        except json.JSONDecodeError:
            return {"answer": raw, "evidence": None, "confidence": "low"}


# ── extract_information ───────────────────────────────────────────────────────

def extract_information(
    papers_file: str,
    question:    str,
    file_dir:    str = "files",
    progress_cb       = None,
) -> dict:
    """
    For each paper in papers_file, ask the LLM the given question and
    save a structured CSV of answers.

    The LLM model is read from LLM_MODEL_LITERATURE env var.

    Args:
        papers_file: path to JSON produced by search_literature().
        question:    free-text question formulated by the mother LLM.
                     e.g. "What proteins are associated with this disease?"
                          "What compounds were tested and what were their IC50 values?"
        file_dir:    used to resolve relative paths; output saved next to papers_file.
        progress_cb: optional callable(str) for progress messages.

    Returns:
        {"n_papers": int, "n_answered": int, "csv_file": str}
    """
    def _p(msg):
        if progress_cb: progress_cb(msg)

    with open(papers_file) as f:
        papers = json.load(f)

    client, client_type, model_name = _get_llm()
    total = len(papers)
    rows  = []

    for i, paper in enumerate(papers):
        _p(f"extracting: {i + 1}/{total}")

        if not paper.get("abstract"):
            extraction = {"answer": None, "evidence": None, "confidence": None}
        else:
            try:
                prompt     = _EXTRACTION_PROMPT.format(
                    question = question,
                    title    = paper.get("title", ""),
                    abstract = paper["abstract"],
                )
                extraction = _call_llm(client, client_type, model_name, prompt)
            except Exception as e:
                extraction = {"answer": None, "evidence": None, "confidence": None, "error": str(e)}

        rows.append({
            "pmid":       paper.get("pmid", ""),
            "title":      paper.get("title", ""),
            "year":       paper.get("year", ""),
            "url":        paper.get("url", ""),
            "source":     paper.get("source", ""),
            "answer":     extraction.get("answer"),
            "evidence":   extraction.get("evidence"),
            "confidence": extraction.get("confidence"),
        })

    df          = pd.DataFrame(rows)
    csv_file    = papers_file.replace(".json", "_extracted.csv")
    df.to_csv(csv_file, index=False)

    n_answered = df["answer"].notna().sum()
    _p(f"done: {n_answered}/{total} papers with answers")

    return {
        "n_papers":   total,
        "n_answered": int(n_answered),
        "csv_file":   csv_file,
    }
