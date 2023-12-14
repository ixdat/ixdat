import re
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent


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
        for m in search_res:
            if 'f"' in m:
                good_matches.append(m)
                count += 1

    if good_matches:
        print("\n### FILE ###", f.relative_to(THIS_DIR))
        for good_match in good_matches:
            print(good_match)

print(count)
