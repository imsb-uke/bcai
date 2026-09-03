# Stage 1 — build the React/Vite bundle
FROM node:20-alpine AS build
WORKDIR /app

COPY package.json ./
# npm install rather than npm ci — no lockfile committed
RUN npm install

COPY . .
RUN npm run build

# Stage 2 — serve the static bundle with nginx
FROM nginx:alpine
COPY --from=build /app/dist      /usr/share/nginx/html
COPY nginx.conf                  /etc/nginx/conf.d/default.conf
EXPOSE 80
