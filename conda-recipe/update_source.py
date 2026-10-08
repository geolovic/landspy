"""Refresh the recipe's source URL/hash from the published PyPI sdist."""
import json
from pathlib import Path
import re
from urllib.request import urlopen

recipe = Path(__file__).with_name('recipe.yaml')
text = recipe.read_text()
version = re.search(r'^  version: "([^"]+)"$', text, re.M).group(1)
with urlopen(f'https://pypi.org/pypi/landspy/{version}/json', timeout=30) as response:
    metadata = json.load(response)
sdist, = [item for item in metadata['urls'] if item['packagetype'] == 'sdist']
text = re.sub(r'^  url: .*$', '  url: ' + sdist['url'], text, flags=re.M)
text = re.sub(r'^  sha256: .*$', '  sha256: ' + sdist['digests']['sha256'], text, flags=re.M)
recipe.write_text(text)
print(f'Updated source for landspy {version} from PyPI.')
