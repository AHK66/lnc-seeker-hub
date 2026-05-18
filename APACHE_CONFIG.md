# Apache 2 Configuration for lnc-seeker-hub

This guide explains how to configure an Apache 2 server as a reverse proxy for the `lnc-seeker-hub` Bokeh application. This setup allows you to access the dashboard via a standard web URL (e.g., `http://your-domain.com/lnc_seeker_server`) and integrates it into your existing web infrastructure.

## 1. Prerequisites

Ensure you have Apache 2 installed and the necessary proxy modules enabled.

### Enable Apache Modules (Linux)
Run the following commands to enable the required modules:
```bash
sudo a2enmod proxy
sudo a2enmod proxy_http
sudo a2enmod proxy_wstunnel
sudo systemctl restart apache2
```

### Enable Apache Modules (Windows)
Uncomment the following lines in your `httpd.conf`:
```apache
LoadModule proxy_module modules/mod_proxy.so
LoadModule proxy_http_module modules/mod_proxy_http.so
LoadModule proxy_wstunnel_module modules/mod_proxy_wstunnel.so
```

## 2. Launching the Bokeh Server

To allow Apache to tunnel requests to Bokeh, you must start the Bokeh server with the `--allow-websocket-origin` flag set to your domain name.

### Production Command
```bash
python -m bokeh serve lnc_seeker_server.py \
    --port 5006 \
    --allow-websocket-origin=your-domain.com \
    --session-token-expiration 3600
```
*Replace `your-domain.com` with your actual domain or IP address.*

## 3. Apache Site Configuration

Add the following configuration to your Apache virtual host file (e.g., `/etc/apache2/sites-available/000-default.conf` on Linux or `extra/httpd-vhosts.conf` on Windows).

```apache
<VirtualHost *:80>
    ServerName your-domain.com

    # Preserve the original Host header
    ProxyPreserveHost On

    # 1. Handle Bokeh WebSockets (Crucial for interactivity)
    # The order matters: WebSocket rules must come before general HTTP rules
    ProxyPass /lnc_seeker_server/ws ws://localhost:5006/lnc_seeker_server/ws
    ProxyPassReverse /lnc_seeker_server/ws ws://localhost:5006/lnc_seeker_server/ws

    # 2. Handle Bokeh HTTP requests
    ProxyPass /lnc_seeker_server http://localhost:5006/lnc_seeker_server
    ProxyPassReverse /lnc_seeker_server http://localhost:5006/lnc_seeker_server

    # (Optional) Redirect the root URL to the application
    # RedirectMatch ^/$ /lnc_seeker_server
</VirtualHost>
```

## 4. How to link from a standard web page

Once Apache is configured, you can link to the lnc-seeker-hub from any other page using a standard anchor tag:

```html
<a href="http://your-domain.com/lnc_seeker_server">Open lnc-seeker Hub</a>
```

## 5. Security and HTTPS (Recommended)

If you are using SSL/TLS (HTTPS), you must change the protocol in the WebSocket proxy from `ws` to `wss`:

```apache
# For HTTPS setup
ProxyPass /lnc_seeker_server/ws wss://localhost:5006/lnc_seeker_server/ws
ProxyPassReverse /lnc_seeker_server/ws wss://localhost:5006/lnc_seeker_server/ws
```

## 6. Keeping the Server Running

For production, it is recommended to use a process manager like `systemd` (Linux) or `nssm` (Windows) to ensure the Bokeh server starts automatically and restarts if it crashes.

### Example Systemd Service (`/etc/systemd/system/lnc-seeker.service`)
```ini
[Unit]
Description=lnc-seeker Bokeh Server
After=network.target

[Service]
Type=simple
User=your-user
WorkingDirectory=/path/to/lnc-seeker-hub
ExecStart=/path/to/lnc-seeker-hub/.venv/bin/python -m bokeh serve lnc_seeker_server.py --port 5006 --allow-websocket-origin=your-domain.com
Restart=always

[Install]
WantedBy=multi-user.target
```
