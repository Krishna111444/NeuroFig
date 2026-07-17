"""
licensing.py — offline-verifiable license keys for NeuroFig.

A license key is a compact, signed token (mini-JWT style):

    base64url(payload) . base64url(HMAC-SHA256(secret, payload))

Only the seller knows the secret, so keys cannot be forged, yet the running app
can verify them **offline** with no database and no per-check network call. The
payload carries the buyer's email, an expiry, and a plan tag.

Seller workflow (after a customer pays):
    NEUROFIG_LICENSE_SECRET=... python licensing.py buyer@lab.edu 365
    -> prints a key; email it to the buyer, who pastes it into the app.

This module has no third-party dependencies and no UI — it is unit-testable.
"""
from __future__ import annotations
import base64
import hashlib
import hmac
import json
import time


def _b64e(b: bytes) -> str:
    return base64.urlsafe_b64encode(b).decode().rstrip("=")


def _b64d(s: str) -> bytes:
    return base64.urlsafe_b64decode(s + "=" * (-len(s) % 4))


def _sign(body: str, secret: str) -> str:
    return _b64e(hmac.new(secret.encode(), body.encode(), hashlib.sha256).digest())


def make_license(email: str, days_valid: int, secret: str,
                 plan: str = "pro", issued_at: int | None = None) -> str:
    """Issue a signed key valid for `days_valid` days from `issued_at` (epoch s)."""
    if not secret:
        raise ValueError("A non-empty secret is required to issue licenses.")
    issued = int(issued_at if issued_at is not None else time.time())
    payload = {"e": email, "iat": issued, "x": issued + days_valid * 86400, "p": plan}
    body = _b64e(json.dumps(payload, separators=(",", ":"), sort_keys=True).encode())
    return f"{body}.{_sign(body, secret)}"


def verify_license(key: str, secret: str, now: int | None = None) -> tuple[bool, str, dict]:
    """Return (valid, reason, info).

    Checks the signature in constant time, then the expiry. `now` (epoch s) is
    injectable for testing. Never raises on bad input — returns a reason instead.
    """
    now = int(now if now is not None else time.time())
    if not key or not secret:
        return False, "missing key or secret", {}
    parts = key.strip().split(".")
    if len(parts) != 2:
        return False, "malformed key", {}
    body, sig = parts
    if not hmac.compare_digest(sig, _sign(body, secret)):
        return False, "invalid signature", {}
    try:
        info = json.loads(_b64d(body))
    except Exception:
        return False, "corrupt payload", {}
    if now > int(info.get("x", 0)):
        return False, "expired", info
    return True, "ok", info


if __name__ == "__main__":  # tiny seller-side key issuer
    import os
    import sys

    secret = os.environ.get("NEUROFIG_LICENSE_SECRET", "")
    if len(sys.argv) < 3 or not secret:
        sys.exit("usage: NEUROFIG_LICENSE_SECRET=... python licensing.py <email> <days> [plan]\n"
                 "  (set the secret in your environment; never commit it)")
    email, days = sys.argv[1], int(sys.argv[2])
    plan = sys.argv[3] if len(sys.argv) > 3 else "pro"
    key = make_license(email, days, secret, plan=plan)
    ok, reason, info = verify_license(key, secret)
    print(key)
    print(f"# {email} · {plan} · valid {days} days · self-check: {reason}")
