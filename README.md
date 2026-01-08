# oligostore
Django, Docker web-server for managing oligonucleotides, such as primers and plasmids

# Local HTTPS Setup (OpenSSL + Docker + NGINX + Django)

This document describes how to enable **HTTPS locally** using **OpenSSL**, with
**NGINX terminating TLS** and **Django running behind Gunicorn** in Docker.

This setup is suitable for:
- local development
- internal staging
- production-like testing

> ⚠️ Browsers will warn about the certificate unless it is explicitly trusted.
> This is expected for self-signed certificates.

---

## Architecture Overview

Browser
└── HTTPS (443)
└── NGINX (TLS termination)
└── HTTP (8000)
└── Gunicorn
└── Django


- TLS is handled **only by NGINX**
- Django and Gunicorn do **not** manage certificates
- Docker internal traffic remains HTTP

---

## 1. Generate Self-Signed Certificates (OpenSSL)

From the project root:

```bash
mkdir -p nginx/certs

openssl req -x509 -nodes -days 365 \
  -newkey rsa:2048 \
  -keyout nginx/certs/localhost.key \
  -out nginx/certs/localhost.crt \
  -subj "/CN=localhost"
```  

### Directory structure:
nginx/certs/
├── localhost.crt   # certificate
└── localhost.key     # private key

# 2. docker-compose nginx service:

nginx:
  image: nginx:latest
  ports:
    - "80:80"
    - "443:443"
  volumes:
    - ./nginx/nginx.conf:/etc/nginx/nginx.conf:ro
    - ./nginx/certs:/etc/nginx/certs:ro
  depends_on:
    - web

# 3. Nginx configuration:

events { }

http {
    include /etc/nginx/mime.types;

    upstream django_app {
        server web:8000;
    }

    # HTTP → HTTPS
    server {
        listen 80;
        server_name localhost;
        return 301 https://$host$request_uri;
    }

    # HTTPS
    server {
        listen 443 ssl;
        server_name localhost;

        ssl_certificate     /etc/nginx/certs/fullchain.pem;
        ssl_certificate_key /etc/nginx/certs/privkey.pem;

        location /static/ {
            alias /app/staticfiles/;
        }

        location /media/ {
            alias /app/media/;
        }

        location / {
            proxy_pass http://django_app;
            proxy_set_header Host $host;
            proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header X-Forwarded-Proto https;
            proxy_set_header X-Forwarded-Port 443;
        }
    }
}

# 4. Django Production Settings
## In settings.py
```bash
SECURE_PROXY_SSL_HEADER = ("HTTP_X_FORWARDED_PROTO", "https")
SECURE_SSL_REDIRECT = True
SESSION_COOKIE_SECURE = True
CSRF_COOKIE_SECURE = True
```
## In .env.production
```bash
DJANGO_DEBUG=False
```
# 5. Start the stack
```bash
docker compose down
docker compose build --no-cache
docker compose up
```

# 6. Verify HTTPS
```bash
openssl s_client -connect localhost:443 -servername localhost
curl.exe -vk https://localhost
```

# 7. Trust the certificate (windows)
```bash
openssl x509 -in nginx/certs/fullchain.pem -out localhost.crt
```