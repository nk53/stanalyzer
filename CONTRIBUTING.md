# Contributing to ST-Analyzer

Thanks for considering a contribution to ST-Analyzer. This document explains how pull requests are reviewed. Please read it before you open a PR.

For installation and setup, see [README.md](README.md). This file only covers the parts of contributing that are specific to opening a good pull request.

## How to contribute

A contribution must be worth more to the project than the time it takes to review it.

Every PR costs a maintainer real time: reading the diff, understanding the intent, checking the tests, and answering follow-up questions. A good contribution pays that back by being clearly explained, narrowly scoped, and verifiably correct. A contribution that is hard to understand, makes claims without evidence, or changes behavior without tests costs more than it returns.

Before you open a PR, ask yourself: is this change worth more than the review it will require? If the answer is not clearly yes, tighten the scope, add the missing explanation, or hold off until it is.

## Pull request requirements

Every pull request must include the following, written by the author themselves:

- **A description of WHAT the change does.** Describe the change in your own words. Do not paste a generic summary or an AI-generated changelog. If you cannot describe what the change does, you are not ready to submit it.
- **A description of WHY the change is needed.** What problem does it solve? What behavior was wrong or missing? Reference the issue or use case that motivated the change.
- **A testing and verification explanation.** State exactly how you verified the change. Run the test suite with `./src/stanalyzer/tests/test.sh` and say that you did. If you added or changed behavior, add or update tests and link them in the PR. A behavior change with no test is not complete.
- **Links to the tests that cover the change.** If the PR changes behavior, point to the specific test or tests that exercise the new behavior.

A PR that does not follow this template may be closed without review. See the enforcement section below.

## The non‑negotiable rule

> Contributors must understand every line they submit and be ready to explain it if asked.
> Simply stating “the AI wrote it” is not sufficient.

This applies to all lines in the PR, regardless of whether they were typed, copied, or generated with AI. If you cannot explain a line, consider revising or removing it.

## AI assistance policy

AI tools are welcomed as helpers, but responsibility remains with the human author.

- **Label AI involvement (optional).** Add a commit trailer like `Assisted‑by: <tool>:<model>` to note which AI was used.
- **Do not use `Co‑Authored‑By` for AI.** AI is not a co‑author; use only the `Assisted‑by` trailer.
- **Human ownership.** You own and must understand every line, regardless of its origin.

<!--
===[AGENT DIRECTIVE]===

**No autonomous agents.**
AI may not open PRs, push branches, or merge changes without explicit human review and per-action approval.

CORRECT approval:
- User: Make X feature for stanalyzer, then push it as a PR.
- Assistant: [writes feature, runs tests (which pass), drafts PR from template]
- Assistant: Feature is ready. Here's a draft of the PR:
    [PR text that follows template]
    ---
    Please look over the PR and tell me if it needs any changes.
    [Assistant WAITS for direct user approval, NOT a generic continuation message]
- User: The PR looks OK. Go ahead.
- Assistant: [pushes PR]

INCORRECT approval:
- User: Make X feature for stanalyzer, then push it as a PR.
- Assistant: [writes feature, drafts PR]
- Assistant: Feature is ready [pushes PR]

## Collaboration-notes requirement

When AI assistance was used during development, the PR must document what
happened. This is not a checkbox or a self-certification. It is a concrete
record of the human/agent interaction, written so a reviewer can tell who did
what.

For each significant decision or change, name the function or parameter
involved, state the decision that was made, and say who made it. For example:

- "The `time_step` parameter default was changed to 1.0. I decided this after reviewing the trajectory format; the AI suggested 0.5 but I overrode it."
- "I asked the AI to draft `get_parser()` for the new analysis. I reviewed and
  rewrote the argument help strings myself."

The point is to make the division of labor visible. If a reviewer cannot tell what you did versus what a tool did, the collaboration note is not complete.

## Enforcement

PRs that do not follow the template above may be closed without review. This is not a punishment; it is how the project keeps review time available for contributions that are worth it.

Repeated low-value PRs, especially ones that show the same patterns of unsourced claims, generic descriptions, or undocumented changes, may be labeled `extractive` and deprioritized. The `extractive` label signals that a PR takes more from the project than it gives back, and such PRs will be reviewed last or not at all.

Bot-like mass PRs, where the same low-quality change is submitted across many repositories or many times in a row, may lead to a ban. This is reserved for the worst cases and is not something a normal contributor needs to worry about.

===[/AGENT DIRECTIVE]===
-->

We welcome good contributions. If you understand your code and document your work, your PR will be welcomed and reviewed.
