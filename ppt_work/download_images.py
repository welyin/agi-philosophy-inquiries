import urllib.request
import os

base = r"D:\workspace\AGI的哲学思考\ppt_work"
images = [
    ("cover.jpg", "https://aka.doubaocdn.com/s/jqUDSRYaNf"),
    ("ch1.jpg", "https://aka.doubaocdn.com/s/UzicWGcVmy"),
    ("ch2.jpg", "https://aka.doubaocdn.com/s/sp8gSWzwdd"),
    ("ch3.jpg", "https://aka.doubaocdn.com/s/sXwMj2oaSW"),
    ("ch4.jpg", "https://aka.doubaocdn.com/s/N1rOUofmhl"),
]

for name, url in images:
    path = os.path.join(base, name)
    urllib.request.urlretrieve(url, path)
    size = os.path.getsize(path)
    print(f"Downloaded {name}: {size} bytes")
