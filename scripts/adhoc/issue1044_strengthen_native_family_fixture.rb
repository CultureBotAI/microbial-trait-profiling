# Purpose: retain the literal fixture edit used for issue #1044, 2026-10-25.
# Input/output: tests/test_three_test_native.py, relative to the repository root.
# This changes only synthetic imported probabilities so selected-family BY is
# distinguishable from padding with an excluded trait; it runs no analysis/tests.

path = "tests/test_three_test_native.py"
source = File.read(path)
old = "    records = list(rows(_native_path(sources)))\n    p_values = [float(row[\"study_p_value\"]) for row in records]\n"
replacement = "    records = list(rows(_native_path(sources)))\n    # Isolate correction-family arithmetic with controlled imported probabilities.\n    # The preservation tests above use the unchanged actual producer outputs.\n    p_values = [0.001 * (index + 1) for index in range(len(records))]\n    native_bh = false_discovery_control(p_values, method=\"bh\")\n    for row, p_value, q_value in zip(records, p_values, native_bh, strict=True):\n        row[\"study_p_value\"] = str(p_value)\n        row[\"study_q_value\"] = str(q_value)\n    write_table(_native_path(sources), records)\n"
abort "Fixture edit already applied; no file changed" if source.include?(replacement)
abort "Expected original fixture literal absent; no file changed" unless source.include?(old)
File.write(path, source.sub(old, replacement))
