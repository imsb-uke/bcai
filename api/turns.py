import asyncio
from typing import Any, AsyncGenerator, Awaitable, Callable, Dict, List, Optional


class RunningTurn:
    """
    One in-flight chat turn, running as an independent asyncio.Task that is not
    tied to any particular HTTP request. Any number of SSE connections — the
    request that started it, plus later reattachments after a session switch,
    a logout/login, a dropped connection, or another browser tab — can all
    subscribe to the same turn concurrently. Each subscriber first replays
    everything emitted so far, then gets live events. Closing one subscriber's
    connection (client abort, tab close) never touches the underlying task —
    only an explicit cancel (see TurnRegistry.cancel) stops it.
    """

    def __init__(self, session_id: str, username: str) -> None:
        self.session_id = session_id
        self.username   = username
        self.buffer: List[Dict[str, Any]] = []
        self.subscribers: List["asyncio.Queue[Optional[Dict[str, Any]]]"] = []
        self.done = False
        self.task: Optional[asyncio.Task] = None

    def publish(self, event: Dict[str, Any]) -> None:
        self.buffer.append(event)
        for q in self.subscribers:
            q.put_nowait(event)

    def finish(self) -> None:
        self.done = True
        for q in self.subscribers:
            q.put_nowait(None)   # sentinel: this subscriber's stream ends

    async def subscribe(self) -> AsyncGenerator[Dict[str, Any], None]:
        queue: "asyncio.Queue[Optional[Dict[str, Any]]]" = asyncio.Queue()
        for event in self.buffer:
            queue.put_nowait(event)
        if self.done:
            queue.put_nowait(None)
        else:
            self.subscribers.append(queue)
        try:
            while True:
                event = await queue.get()
                if event is None:
                    return
                yield event
        finally:
            if queue in self.subscribers:
                self.subscribers.remove(queue)


class TurnRegistry:
    """
    Tracks currently-running turns, keyed by session_id — client-generated
    UUIDs (crypto.randomUUID() in ChatPage.tsx), globally unique, so no
    per-username namespacing is needed here.
    """

    def __init__(self) -> None:
        self._turns: Dict[str, RunningTurn] = {}

    def get(self, session_id: str) -> Optional[RunningTurn]:
        turn = self._turns.get(session_id)
        return turn if turn is not None and not turn.done else None

    def is_running(self, session_id: str) -> bool:
        return self.get(session_id) is not None

    def start(
        self,
        session_id: str,
        username: str,
        run: Callable[[RunningTurn], Awaitable[None]],
    ) -> RunningTurn:
        turn = RunningTurn(session_id, username)
        self._turns[session_id] = turn

        async def _runner() -> None:
            try:
                await run(turn)
            finally:
                turn.finish()
                if self._turns.get(session_id) is turn:
                    del self._turns[session_id]

        turn.task = asyncio.create_task(_runner())
        return turn

    def cancel(self, session_id: str) -> bool:
        turn = self.get(session_id)
        if turn is None or turn.task is None:
            return False
        turn.task.cancel()
        return True
