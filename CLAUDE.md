# CLAUDE.md

## Commit Rule: Every Step

After each completed step, commit and push:

```bash
git add -A
git commit -m "Prefix: description"
git pull --rebase
git push
```

Prefixes: `Add:` / `Update:` / `Complete:` / `Fix:` / `Archive:`
Always include: `Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>`
