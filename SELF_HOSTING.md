# Self-hosting NeuroFig on your own server

A step-by-step runbook to run NeuroFig on a personal Ubuntu VPS behind a real
domain with HTTPS, surviving reboots and crashes. Everything here is copy-paste.

Estimated cost: a small VPS (2 GB RAM) runs ~$5–12/month (Hetzner, DigitalOcean,
Linode, Vultr). 1 GB works but 2 GB is comfortable for matplotlib + SciPy.

Legend: run commands as your normal sudo user unless a step says otherwise.
Replace `neurofig.example.com` with your real domain and `<...>` placeholders.

---

## 0. What you need first
- An Ubuntu 22.04 or 24.04 server with SSH access (a "VPS" or "droplet").
- A domain name you control (for HTTPS). Optional but strongly recommended.
- Your payment link (Lemon Squeezy / Razorpay) — only needed when you switch on
  payments; you can deploy free-tier first and add it later.

---

## 1. Point your domain at the server
In your domain registrar's DNS, add an **A record**:

| Type | Name | Value |
|------|------|-------|
| A | neurofig (or @) | your server's IP address |

DNS can take a few minutes to propagate. Check with `ping neurofig.example.com`.

---

## 2. Basic server setup
```bash
sudo apt update && sudo apt -y upgrade
sudo apt -y install python3-venv python3-pip git nginx openbabel
```
`openbabel` here installs the `obabel` command system-wide — NeuroFig auto-detects
it (via `shutil.which`), so the PDB→PDBQT tab just works, no conda needed.

Create a dedicated, unprivileged user to run the app (good security hygiene):
```bash
sudo adduser --system --group --home /opt/neurofig neurofig
```

---

## 3. Get the code onto the server
If your code is on GitHub:
```bash
sudo -u neurofig git clone https://github.com/<you>/neurofig.git /opt/neurofig/app
```
No GitHub? Copy it up from your machine instead:
```bash
# run this on YOUR computer, from the project folder
rsync -av --exclude '.venv' --exclude '__pycache__' ./ \
  <you>@<server-ip>:/tmp/neurofig/
# then on the server:
sudo mkdir -p /opt/neurofig/app && sudo cp -r /tmp/neurofig/* /opt/neurofig/app/
sudo chown -R neurofig:neurofig /opt/neurofig
```

---

## 4. Create the Python environment
```bash
sudo -u neurofig python3 -m venv /opt/neurofig/venv
sudo -u neurofig /opt/neurofig/venv/bin/pip install --upgrade pip
sudo -u neurofig /opt/neurofig/venv/bin/pip install -r /opt/neurofig/app/requirements.txt
```
Verify the engine works before going further:
```bash
cd /opt/neurofig/app && sudo -u neurofig /opt/neurofig/venv/bin/python -m tests.test_core
```
You should see the tests pass. (Open Babel test: `python -m tests.test_pdbqt`.)

---

## 5. Configure secrets
Create `/opt/neurofig/app/.streamlit/secrets.toml` on the server only (it is
gitignored and must never be committed):
```bash
sudo -u neurofig mkdir -p /opt/neurofig/app/.streamlit
sudo -u neurofig tee /opt/neurofig/app/.streamlit/secrets.toml >/dev/null <<'EOF'
# Long random string, used to sign AND verify license keys. Keep it secret.
NEUROFIG_LICENSE_SECRET = "REPLACE_ME"
# Your checkout link. Leave "" to hide the Buy button until you're ready to sell.
NEUROFIG_BUY_URL = ""
EOF
sudo chmod 600 /opt/neurofig/app/.streamlit/secrets.toml
```
Generate the secret once and paste it in:
```bash
python3 -c "import secrets; print(secrets.token_urlsafe(32))"
```
Also make the app's server config headless-friendly:
```bash
sudo -u neurofig tee /opt/neurofig/app/.streamlit/config.toml >/dev/null <<'EOF'
[server]
headless = true
address = "127.0.0.1"
port = 8501
maxUploadSize = 50
enableXsrfProtection = true
gatherUsageStats = false
[browser]
gatherUsageStats = false
EOF
```

---

## 6. Run it as a service (survives reboots and crashes)
Create the systemd unit:
```bash
sudo tee /etc/systemd/system/neurofig.service >/dev/null <<'EOF'
[Unit]
Description=NeuroFig Streamlit app
After=network.target

[Service]
Type=simple
User=neurofig
Group=neurofig
WorkingDirectory=/opt/neurofig/app
Environment=MPLBACKEND=Agg
ExecStart=/opt/neurofig/venv/bin/streamlit run app.py
Restart=always
RestartSec=3

[Install]
WantedBy=multi-user.target
EOF

sudo systemctl daemon-reload
sudo systemctl enable --now neurofig
sudo systemctl status neurofig --no-pager
```
`Restart=always` means it comes back automatically if it crashes or the server
reboots. Logs: `journalctl -u neurofig -f`.

At this point the app is running on `127.0.0.1:8501` — reachable only from the
server itself. nginx will expose it to the world next.

---

## 7. Put nginx in front (HTTP)
Streamlit uses WebSockets, so the proxy config below includes the Upgrade headers
it needs — miss those and the page loads but never connects.
```bash
sudo tee /etc/nginx/sites-available/neurofig >/dev/null <<'EOF'
server {
    listen 80;
    server_name neurofig.example.com;

    client_max_body_size 50M;

    location / {
        proxy_pass http://127.0.0.1:8501;
        proxy_http_version 1.1;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        # WebSocket support (required by Streamlit)
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_read_timeout 86400;
    }
}
EOF

sudo ln -sf /etc/nginx/sites-available/neurofig /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default
sudo nginx -t && sudo systemctl reload nginx
```
Now visit `http://neurofig.example.com` — you should see NeuroFig.

---

## 8. Turn on HTTPS (free, automatic)
```bash
sudo apt -y install certbot python3-certbot-nginx
sudo certbot --nginx -d neurofig.example.com
```
Certbot edits your nginx config to serve HTTPS and sets up auto-renewal. Choose
"redirect HTTP to HTTPS" when asked. Your site is now `https://...`.

---

## 9. Lock down the firewall
```bash
sudo ufw allow OpenSSH
sudo ufw allow 'Nginx Full'
sudo ufw --force enable
```
Port 8501 is never exposed — only 22 (SSH), 80, and 443 are open.

---

## 10. Selling: issue a license key after each payment
Free tier (watermarked previews) works with `NEUROFIG_BUY_URL = ""`. When you're
ready to charge:
1. Set `NEUROFIG_BUY_URL` in `secrets.toml` to your checkout link, then
   `sudo systemctl restart neurofig`.
2. After a customer pays, mint them a key on the server:
```bash
cd /opt/neurofig/app
sudo -u neurofig env NEUROFIG_LICENSE_SECRET="<same secret as secrets.toml>" \
  /opt/neurofig/venv/bin/python licensing.py buyer@lab.edu 365
```
   This prints a signed key valid 365 days — email it to the buyer, who pastes it
   into the sidebar to unlock clean vector exports. Verification is offline; no
   database needed. (Automate later with a payment webhook — manual is fine to start.)

---

## 11. Updating the app later
```bash
cd /opt/neurofig/app
sudo -u neurofig git pull                       # or rsync the new files up
sudo -u neurofig /opt/neurofig/venv/bin/pip install -r requirements.txt
sudo systemctl restart neurofig
```

---

## 12. Keep the promise you make in the app
The sidebar tells users their data is not stored after the session and not used
for training. Honour it:
- NeuroFig already holds uploads only in memory for the request — it writes no
  user data to disk (the PDBQT tool uses an auto-deleted temp dir).
- Don't add logging that records uploaded file contents.
- Keep the server patched: `sudo apt update && sudo apt upgrade` periodically.

---

## Troubleshooting
| Symptom | Fix |
|---|---|
| Page loads but spins forever | WebSocket headers missing in nginx (step 7) — check `Upgrade`/`Connection` lines, `sudo nginx -t`, reload |
| 502 Bad Gateway | app not running: `sudo systemctl status neurofig`, `journalctl -u neurofig -e` |
| PDB→PDBQT says Open Babel missing | `which obabel` should print `/usr/bin/obabel`; if not, `sudo apt install openbabel`. Or set `Environment=NEUROFIG_OBABEL=/usr/bin/obabel` in the service file |
| "License key rejected" for a key you issued | the secret used to *issue* must equal `NEUROFIG_LICENSE_SECRET` in `secrets.toml` |
| Out-of-memory / killed | use a 2 GB server, or add swap: `sudo fallocate -l 2G /swapfile && sudo chmod 600 /swapfile && sudo mkswap /swapfile && sudo swapon /swapfile` |
| Uploads rejected as too large | raise `maxUploadSize` (config.toml) and `client_max_body_size` (nginx) together |

---

## One-glance architecture
```
Browser ──HTTPS──> nginx (:443)  ──proxy──>  Streamlit (127.0.0.1:8501)
                     │                          │
                 Let's Encrypt              NeuroFig engine
                 certificate                 (calls system `obabel` for PDBQT)
```
That's the whole stack: one small server, one service, one reverse proxy.
```
