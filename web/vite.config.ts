import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

export default defineConfig({
  plugins: [react()],
  server: {
    // Dev proxy mirrors what nginx does in prod — /api/* → backend on 8000
    proxy: {
      '/api': {
        target      : 'http://localhost:8000',
        rewrite     : (path) => path.replace(/^\/api/, ''),
        changeOrigin: true,
      },
    },
  },
})
