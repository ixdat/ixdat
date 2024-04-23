"""This file will look for f-strings in exceptions messages of any .py files **in the folder it is
placed in** or its sub folders."""

import re
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent


# Look for the 'raise' keyword and capture anything in parentheses (non greedily) after that keyword
RAISE = re.compile("raise .*?\([\s\S]*?\)")


count = 0
for f in THIS_DIR.glob("**/*"):
    if f.suffix != ".py":
        continue

    with open(f) as file_:
        content = file_.read()

    search_res = RAISE.findall(content)
    good_matches = []
    if search_res:
        # Look for the f-string marker in the found exceptions strings
        for m in search_res:
            if 'f"' in m:
                good_matches.append(m)
                count += 1

    if good_matches:
        print("\n### FILE ###", f.relative_to(THIS_DIR))
        for good_match in good_matches:
            print(good_match)

print(count)
