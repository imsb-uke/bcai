import { useAuth }  from './store';
import AuthPage    from './components/AuthPage';
import ChatPage    from './components/ChatPage';

export default function App() {
  const { token } = useAuth();
  // Token presence is the only routing decision — no router library needed
  return token ? <ChatPage /> : <AuthPage />;
}
