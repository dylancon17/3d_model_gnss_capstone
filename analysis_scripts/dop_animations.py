import os
import json
import argparse

IMAGE_EXTS = (".jpg", ".jpeg", ".png", ".gif", ".webp")

parser = argparse.ArgumentParser(description="Generate HTML timelapse from images")
parser.add_argument("folder", help="Path to folder containing images")
args = parser.parse_args()

IMAGE_FOLDER = args.folder.rstrip("/\\")

frames = []

#Load images
if not os.path.isdir(IMAGE_FOLDER):
    raise FileNotFoundError(f"Folder not found: {IMAGE_FOLDER}")

OUTPUT_HTML = IMAGE_FOLDER + "\\timelapse.html"

images = sorted(
    f for f in os.listdir(IMAGE_FOLDER)
    if f.lower().endswith(IMAGE_EXTS)
)

for img in images:
    frames.append(f"{IMAGE_FOLDER}/{img}")

if not frames:
    raise RuntimeError("No images found in the specified folder")

#Generate HTML
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Image Timelapse</title>
<style>
    body {{
        background: #111;
        color: #eee;
        font-family: sans-serif;
        text-align: center;
    }}
    img {{
        max-width: 90%;
        max-height: 70vh;
        margin-top: 20px;
    }}
    input[type=range] {{
        width: 80%;
    }}
    button {{
        padding: 8px 16px;
        margin: 10px;
        font-size: 16px;
    }}
    .label {{
        margin-top: 10px;
    }}
</style>
</head>
<body>

<h2>Timelapse Viewer</h2>

<img id="frame" src="{frames[0]}">

<div>
    <input type="range" id="slider" min="0" max="{len(frames)-1}" value="0">
</div>

<div class="label">
    Speed: <span id="speedValue">5.0</span> sec / frame
</div>
<div>
    <input type="range" id="speedSlider" min="0.1" max="10" step="0.1" value="5">
</div>

<div>
    <button onclick="play()">▶ Play</button>
    <button onclick="pause()">⏸ Pause</button>
</div>

<script>
const frames = {json.dumps(frames)};
const img = document.getElementById("frame");
const slider = document.getElementById("slider");
const speedSlider = document.getElementById("speedSlider");
const speedValue = document.getElementById("speedValue");

let index = 0;
let timer = null;
let frameDelay = 5000;

function showFrame(i) {{
    index = i;
    img.src = frames[index];
    slider.value = index;
}}

slider.addEventListener("input", () => {{
    showFrame(parseInt(slider.value));
}});

speedSlider.addEventListener("input", () => {{
    speedValue.textContent = speedSlider.value;
    frameDelay = parseFloat(speedSlider.value) * 1000;

    if (timer) {{
        pause();
        play();
    }}
}});

function step() {{
    index = (index + 1) % frames.length;
    showFrame(index);
}}

function play() {{
    if (timer) return;
    timer = setInterval(step, frameDelay);
}}

function pause() {{
    clearInterval(timer);
    timer = null;
}}
</script>

</body>
</html>
"""

#Write to html file
with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
    f.write(html)

print(f"Created {OUTPUT_HTML} with {len(frames)} frames from '{IMAGE_FOLDER}'")