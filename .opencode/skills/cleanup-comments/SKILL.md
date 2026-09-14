---
name: cleanup-comments
description: Clean up code comments to follow minimal commenting guidelines. Use when the user says "clean up comments", "fix comments", "comments are too verbose", "too many comments", or wants to enforce a concise commenting style. Enforces: 1-line max, explain WHY not WHAT, keywords over prose, URLs only if they'll survive 20 years.
---

# Cleanup Comments

Minimal commenting style. Comments document *why* a line exists, not *what* it does.

## Rules

1. **One line.** If you think you need 2+ lines, you don't.
2. **Must match a real line.** Never comment a deletion—put that in the commit message.
3. **WHY, not WHAT.** The line's purpose is usually obvious. The reason it's *needed* usually isn't.
4. **Keywords > prose.** `// HACK: legacy API` > `// This is a hack because we're using a legacy API that doesn't support the new format`.
5. **URLs only if they'll survive 20 years.** Official docs: ok. StackOverflow: no. Random blog: no.
6. **Don't over-explain.** The comment only needs to point a future reader in the right direction, not hand-hold them through the full context.

## What to do

1. Read the file.
2. For each comment:
   - If it explains WHAT the line does (redundant) → **delete**.
   - If it's 2+ lines and can be condensed to 1 → **condense**.
   - If it's prose and keywords would do → **convert to keywords**.
   - If it references a URL that will rot → **replace with project-internal context** (ticket #, bug ID, function name).
   - If it comments a deleted line → **delete**.
   - If it documents a non-obvious WHY and is already 1 line → **keep**.
3. After cleanup, every remaining comment must:
   - Be exactly 1 line.
   - Correspond to an actual line in the file.
   - Explain why the line exists (not what it does).
