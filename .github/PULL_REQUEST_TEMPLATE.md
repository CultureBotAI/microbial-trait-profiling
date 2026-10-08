## Summary

<!-- What does this PR do and why? -->

## Issue

<!-- Link the issue this PR addresses: "Closes #N" or "Refs #N" -->

## Checklist

- [ ] `just check` passes
- [ ] New code has tests
- [ ] No `noqa`, `type: ignore`, or `nosec` comments added
- [ ] Result claims follow [the interpretation guide](../methods/biological_interpretation.md): annotation units, direction/background, effect, uncertainty, support, inference level, scope and provenance are explicit (or no result claims changed)
- [ ] Captions, legends and generated narratives were reviewed; affected claims and evidence/corrections are recorded in [the audit](../methods/biological_interpretation_audit.md)
- [ ] Mechanisms remain hypotheses unless supported; zero signal and shared-data agreement are not presented as absence or independent validation
- [ ] Agent guidance was regenerated with `just sync-agent-docs` if `AGENTS.md` changed
