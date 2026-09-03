import { StrictMode }          from 'react';
import { createRoot }          from 'react-dom/client';
import { AuthProvider, applyTheme, getTheme } from './store';
import App                     from './App';
import './styles/global.css';

// Apply saved theme before first paint to avoid flash
applyTheme(getTheme());

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    <AuthProvider>
      <App />
    </AuthProvider>
  </StrictMode>,
);
