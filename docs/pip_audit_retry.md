# pip-audit retry wrapper

`scripts/run_pip_audit.py` is the serial wrapper used by `just security`,
`just security-linux` and the pip-audit pre-commit hook. These commands retain
their audit scope and arguments; the wrapper does not ignore advisories or
replace an incomplete audit with a successful result.

## Exit and retry policy

A completed successful pip-audit invocation exits 0. Known vulnerabilities
retain pip-audit's exit 1, and argument errors retain argparse's exit 2. The
wrapper exits **3** after exhausting transient failures or refusing a recognized
network error without retry; this distinguishes an incomplete network audit
from findings and invalid arguments (#1163, #2073). Other
nonzero pip-audit results are propagated unchanged. String-valued SystemExit
messages are printed once and exit 1.

There are at most four attempts by default, separated by 2, 4 and 8 seconds.
Recognized connection, timeout, incomplete-read and protocol failures, or
HTTP 408/429/5xx responses, may retry. Typed TLS closures `SSLEOFError` and
`SSLZeroReturnError` also retry, whether bare or nested in requests/urllib3
exceptions. Certificate and hostname verification (`SSLCertVerificationError`
and `CertificateError`), proxy and redirect errors do not retry. HTTP errors
without a recognized transient status do not retry. An explicit permanent cause
takes precedence over a broader connection wrapper. An unclassified bare
`ssl.SSLError` is not guessed transient and is refused without retry; a generic
requests connection error retains the connection retry policy unless its typed
chain contains a permanent cause. Diagnostic words do not classify certificates.
Structured urllib3 `reason` exceptions, exception-valued arguments and
cause/context chains are checked with cycle detection.

Recognized network refusals, including non-transient or unclassified HTTP
errors, print `[pip-audit] NETWORK ERROR (not retried):` followed by typed detail
and exit 3 without a traceback. Unknown programming errors and service errors
without a recognized network cause still raise; they are not guessed transient
or relabeled as network failures.

Console stdout and stderr are forwarded during each attempt. Retry failures
retain their output, and existing console handlers retain their formatters.
The wrapper temporarily routes root and pip-audit console handlers, preserving
other handlers; a temporary console handler supplies diagnostics if none is
present. Wrapper-owned routing is restored after each attempt. A copy of stderr
is retained for the attempt's finding diagnostic; stdout is forwarded without
buffering a copy. Stderr text never determines whether an audit succeeded or
whether a network failure merits a retry.

## Adapter decision (#1164, #1185)

The adapter remains in process for this change. It can classify typed transport
causes, including a certificate or proxy error below a generic connection error.
A subprocess would isolate argv/logging state and avoid rebinding private names,
but pip-audit currently exposes no separate network exit code: retry decisions
would depend on changing stderr wording. A missed text match would still fail
closed, but would lose retries. The existing synthetic-query tests exercise the
real installed CLI and pin the typed retry path. This is the reason for retaining
that path, rather than claiming subprocess isolation is unnecessary.

Before an audit attempt, the wrapper checks whether the top-level `pip_audit`
package is present. If it is not installed in this environment, the CLI reports
an installation error advising `just sync`, then exits 1 before private imports.
For an installed package, the wrapper imports the CLI outside stream capture and
checks its callable and exception aliases. Missing or incompatible private APIs
raise a compatibility error naming the installed pip-audit version; the CLI
reports that error and exits 1. The interface module is prepared once per real
run, rather than imported for each attempt. The CLI exception alias and argv
object are restored to their actual saved values, including caller-owned or
nested bindings. Tests compare the installed version to `uv.lock`, check the
private names, and exercise retry/refusal behavior against the real CLI with
synthetic responses. A dependency update must pass those tests; matching symbol
names alone is not evidence that the upstream exception behavior is unchanged.

This wrapper is intended for serial CLI use. Its temporary argv and logging
changes are process-global; it is not an API for concurrent calls from threads.
No live vulnerability-feed result is established by the deterministic test suite.
