# oligostore

`oligostore` is a Django + Docker application for managing oligonucleotides, including primers, primer pairs, projects, and sequence files.

## What the Tool Does

- Manage primers and primer pairs
- Group laboratory work into projects
- Upload and associate sequence files with projects
- Run sequence and primer-binding analysis workflows
- Export selected primers and primer pairs

## Quick Start (Development)

### 1. Create `.env`

Create a `.env` file in the project root with at least:

```env
DJANGO_SECRET_KEY=change-me
DJANGO_DEBUG=True
DJANGO_ALLOWED_HOSTS=localhost,127.0.0.1
POSTGRES_DB=oligostore
POSTGRES_USER=oligostore
POSTGRES_PASSWORD=oligostore
POSTGRES_HOST=db
POSTGRES_PORT=5432
```

### 2. Start services

```bash
docker compose up --build
```

### 3. Run migrations and create an admin user

In a second terminal:

```bash
docker compose exec web python manage.py migrate
docker compose exec web python manage.py createsuperuser
```

### 4. Open the app

Open `http://localhost:8000`.

## How to Use oligostore

After starting the stack and logging in:

1. Create or import primers from `Primer List`.
2. Create primer pairs from existing primers.
3. Create a project from `Projects`.
4. Attach primer pairs and sequence files on the project dashboard.
5. Run sequence analysis or primer-binding analysis.
6. Download selected primers, primer pairs, or project sequence files.

## Local HTTPS Setup (OpenSSL + Docker + NGINX + Django)

Use this for production-like local HTTPS.

### Architecture

```text
Browser
  -> HTTPS (443)
  -> NGINX (TLS termination)
  -> HTTP (8000)
  -> Gunicorn
  -> Django
```

- TLS is terminated by NGINX.
- Django and Gunicorn run behind NGINX over internal HTTP.

### 1. Generate a self-signed certificate

From the project root:

```bash
mkdir -p nginx/certs
openssl req -x509 -nodes -days 365 -newkey rsa:2048 \
  -keyout nginx/certs/localhost.key \
  -out nginx/certs/localhost.crt \
  -subj "/CN=localhost"
```

Expected files:

```text
nginx/certs/localhost.crt
nginx/certs/localhost.key
```

### 2. Prepare production env files

`docker-compose-production.yml` uses `.env.production` (and currently `.env` for `worker`), so ensure both files exist with required Django/PostgreSQL values.

At minimum in `.env.production`:

```env
DJANGO_SECRET_KEY=change-me
DJANGO_DEBUG=False
DJANGO_ALLOWED_HOSTS=localhost,127.0.0.1
DJANGO_CSRF_TRUSTED_ORIGINS=https://localhost
POSTGRES_DB=oligostore
POSTGRES_USER=oligostore
POSTGRES_PASSWORD=oligostore
POSTGRES_HOST=db
POSTGRES_PORT=5432
```

### 3. Start the production stack

```bash
docker compose -f docker-compose-production.yml up --build
```

### 4. Verify HTTPS

```bash
openssl s_client -connect localhost:443 -servername localhost
curl -vk https://localhost
```

### 5. Trust certificate on Windows (optional)

Export and import the certificate into your trusted root store:

```bash
openssl x509 -in nginx/certs/localhost.crt -out localhost.crt
```

If the certificate is not trusted, browser warnings are expected for self-signed certificates.
