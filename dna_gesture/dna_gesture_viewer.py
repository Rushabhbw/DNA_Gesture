import cv2
import mediapipe as mp
from Bio import SeqIO
from Bio.Seq import Seq
import time
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import os

# -------------------- File Selection --------------------
root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename(title="Select a FASTA File", filetypes=[("FASTA files", "*.fasta *.fa *.txt")])
if not file_path:
    print("No file selected. Exiting.")
    exit()

record = next(SeqIO.parse(file_path, "fasta"))
sequence = str(record.seq).upper()
gene_name = record.id
gene_desc = record.description

# -------------------- Utility Functions --------------------
def calculate_gc(seq):
    g = seq.count('G')
    c = seq.count('C')
    return round(((g + c) / len(seq)) * 100, 2)

def calculate_at(seq):
    a = seq.count('A')
    t = seq.count('T')
    return round(((a + t) / len(seq)) * 100, 2)

def generate_gc_plot(sequence, chunk_size):
    gc_values, x_positions = [], []
    for i in range(0, len(sequence), chunk_size):
        chunk = sequence[i:i+chunk_size]
        if len(chunk) < chunk_size: continue
        gc_values.append(calculate_gc(chunk))
        x_positions.append(i)
    plt.figure(figsize=(6, 2))
    plt.plot(x_positions, gc_values, color='green', marker='o')
    plt.xlabel("Position")
    plt.ylabel("GC%")
    plt.ylim(0, 100)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("gc_plot.png")
    plt.close()

def generate_codon_bar_chart(chunk):
    codons = [chunk[i:i+3] for i in range(0, len(chunk)-2, 3)]
    codons = [c for c in codons if len(c) == 3]
    freq = Counter(codons)
    plt.figure(figsize=(4, 3))
    if freq:
        top_codons = freq.most_common(10)
        labels = [f"{codon} ({count})" for codon, count in top_codons]
        counts = [count for codon, count in top_codons]
        plt.barh(labels[::-1], counts[::-1], color='skyblue')
        plt.xlabel("Frequency")
        plt.title("Top 10 Codons")
    else:
        plt.text(0.5, 0.5, "No codons found", ha='center')
    plt.tight_layout()
    plt.savefig("codon_bar.png")
    plt.close()

def find_orfs(dna_chunk):
    start_codon = 'ATG'
    stop_codons = {'TAA', 'TAG', 'TGA'}
    orfs = []
    for frame in range(3):
        seq = dna_chunk[frame:]
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon == start_codon:
                for j in range(i + 3, len(seq) - 2, 3):
                    stop_codon = seq[j:j+3]
                    if stop_codon in stop_codons:
                        orf_seq = seq[i:j+3]
                        if len(orf_seq) >= 30:
                            aa = str(Seq(orf_seq).translate())
                            orfs.append((i + frame, j + 3 + frame, aa))
                        break
    return orfs

def save_summary(chunk, aa_seq, start):
    orfs = find_orfs(chunk)
    with open("gene_summary.txt", "w") as f:
        f.write(f"Gene Name: {gene_name}\n")
        f.write(f"Description: {gene_desc}\n")
        f.write(f"Length: {len(sequence)} bp\n")
        f.write(f"GC%: {calculate_gc(sequence)}\n")
        f.write(f"AT%: {calculate_at(sequence)}\n")
        f.write(f"Current View (position {start}-{start+len(chunk)}):\n")
        f.write(f"DNA: {chunk}\n")
        f.write(f"AA: {aa_seq}\n\n")
        f.write("Detected ORFs:\n")
        if orfs:
            for idx, (orf_start, orf_end, aa) in enumerate(orfs):
                f.write(f"  ORF {idx+1}: {orf_start}-{orf_end} | {aa}\n")
        else:
            f.write("  No ORFs found.\n")
    print("[✔] gene_summary.txt saved with ORF info.")

def export_html_report(chunk, aa_seq, start):
    orfs = find_orfs(chunk)
    orf_section = ""
    if orfs:
        for idx, (orf_start, orf_end, aa) in enumerate(orfs):
            orf_section += f"<li><strong>ORF {idx+1}:</strong> {orf_start}-{orf_end} → {aa}</li>\n"
    else:
        orf_section = "<li>No ORFs found in this chunk.</li>"
    html = f"""
    <html>
    <head><title>Gene Report</title>
    <style>body {{ font-family: Arial; margin: 40px; }}
    h1 {{ color: #2E8B57; }} .mono {{ font-family: monospace; white-space: pre-wrap; word-wrap: break-word; }}</style></head>
    <body>
    <h1>Gene Summary Report</h1>
    <p><strong>Gene Name:</strong> {gene_name}<br>
    <strong>Description:</strong> {gene_desc}<br>
    <strong>Length:</strong> {len(sequence)} bp<br>
    <strong>GC%:</strong> {calculate_gc(sequence)}%<br>
    <strong>AT%:</strong> {calculate_at(sequence)}%</p>
    <h2>Current Chunk (Position {start}-{start+len(chunk)})</h2>
    <div class="mono">{chunk}</div>
    <p><strong>Amino Acid:</strong><br><div class="mono">{aa_seq}</div></p>
    <h2>Detected ORFs</h2><ul>{orf_section}</ul>
    <h2>GC Content Plot</h2><img src="gc_plot.png">
    <h2>Top 10 Codon Usage</h2><img src="codon_bar.png">
    </body></html>
    """
    with open("gene_report.html", "w", encoding="utf-8") as f:
        f.write(html)
    print("[✔] gene_report.html exported with ORFs.")

# -------------------- Setup --------------------
mp_hands = mp.solutions.hands
hands = mp_hands.Hands(max_num_hands=1)
mp_draw = mp.solutions.drawing_utils
tip_ids = [4, 8, 12, 16, 20]
cap = cv2.VideoCapture(0)

start = 0
chunk_size = 60
mode = "Paused"
last_action_time = 0
base_colors = {'A': (0, 102, 255), 'T': (255, 0, 127), 'G': (0, 255, 0), 'C': (153, 51, 255), 'N': (160, 160, 160)}
show_guide = True

generate_gc_plot(sequence, chunk_size)

# -------------------- Main Loop --------------------
while True:
    success, img = cap.read()
    if not success:
        break
    img = cv2.flip(img, 1)
    img = cv2.resize(img, (1280, 960))
    img_rgb = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    results = hands.process(img_rgb)

    fingers = [0, 0, 0, 0, 0]
    current_time = time.time()

    if results.multi_hand_landmarks:
        handLms = results.multi_hand_landmarks[0]
        lmList = [(int(lm.x * img.shape[1]), int(lm.y * img.shape[0])) for lm in handLms.landmark]
        if lmList:
            fingers[0] = int(lmList[4][0] > lmList[3][0])
            for id in range(1, 5):
                fingers[id] = int(lmList[tip_ids[id]][1] < lmList[tip_ids[id] - 2][1])
        mp_draw.draw_landmarks(img, handLms, mp_hands.HAND_CONNECTIONS)

    finger_count = fingers.count(1)
    if current_time - last_action_time > 1:
        if finger_count == 1:
            start += chunk_size
            mode = "Next"
        elif finger_count == 2:
            start = max(0, start - chunk_size)
            mode = "Previous"
        elif finger_count == 3:
            mode = f"GC: {calculate_gc(sequence[start:start+chunk_size])}%"
        elif finger_count == 4:
            export_html_report(sequence[start:start+chunk_size], Seq(sequence[start:start+chunk_size]).translate(to_stop=True), start)
            mode = "Exported HTML"
        elif finger_count == 5:
            save_summary(sequence[start:start+chunk_size], Seq(sequence[start:start+chunk_size]).translate(to_stop=True), start)
            mode = "Saved TXT"
        last_action_time = current_time

    chunk = sequence[start:start + chunk_size]
    aa_seq = str(Seq(chunk).translate(to_stop=True))
    orfs = find_orfs(chunk)

    # Panel box
    panel_x, panel_y = 650, 10
    cv2.rectangle(img, (panel_x, panel_y), (950, 300), (255, 255, 255), -1)
    cv2.putText(img, f"{gene_name}", (panel_x + 5, panel_y + 20), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), 1)
    cv2.putText(img, f"Mode: {mode}", (panel_x + 5, panel_y + 40), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (150,0,150), 1)
    cv2.putText(img, f"Pos: {start}-{start + chunk_size}", (panel_x + 5, panel_y + 60), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (0,100,255), 1)

    y_base = panel_y + 80
    max_chars_per_line = 25
    for row in range((len(chunk) + max_chars_per_line - 1) // max_chars_per_line):
        x = panel_x + 5
        y = y_base + row * 20
        line_seq = chunk[row * max_chars_per_line: (row + 1) * max_chars_per_line]
        for base in line_seq:
            color = base_colors.get(base.upper(), (0, 0, 0))
            cv2.putText(img, base, (x, y), cv2.FONT_HERSHEY_SIMPLEX, 0.5, color, 1)
            x += 14

    y = y_base + 3 * 20 + 10
    cv2.putText(img, f"AA: {aa_seq}", (panel_x + 5, y), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (0,150,0), 1)
    cv2.putText(img, "ORFs:", (panel_x + 5, y + 20), cv2.FONT_HERSHEY_SIMPLEX, 0.4, (0,100,100), 1)
    for i, (orf_start, orf_end, aa) in enumerate(orfs[:2]):
        cv2.putText(img, f"{orf_start}-{orf_end}: {aa}", (panel_x + 5, y + 40 + i*20), cv2.FONT_HERSHEY_SIMPLEX, 0.35, (50,50,50), 1)

    # Gesture Guide
    if show_guide:
        cv2.rectangle(img, (10, 10), (420, 190), (255, 255, 255), -1)
        cv2.putText(img, "Gesture Controls:", (20, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 0, 0), 2)
        cv2.putText(img, "1 finger: Next Chunk", (20, 55), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (30,30,30), 1)
        cv2.putText(img, "2 fingers: Previous Chunk", (20, 75), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (30,30,30), 1)
        cv2.putText(img, "3 fingers: Show GC%", (20, 95), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (30,30,30), 1)
        cv2.putText(img, "4 fingers: Export HTML", (20, 115), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (30,30,30), 1)
        cv2.putText(img, "5 fingers: Save Summary", (20, 135), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (30,30,30), 1)
        cv2.putText(img, "Press 'h' to hide/show guide", (20, 165), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (100, 0, 0), 1)

    # GC and Codon Images
    if os.path.exists("gc_plot.png"):
        gc_img = cv2.imread("gc_plot.png")
        gc_img = cv2.resize(gc_img, (400, 120))
        img[740:860, 10:410] = gc_img
    if os.path.exists("codon_bar.png"):
        codon_img = cv2.imread("codon_bar.png")
        codon_img = cv2.resize(codon_img, (400, 150))
        img[740:890, 420:820] = codon_img
    else:
        generate_codon_bar_chart(chunk)

    cv2.imshow("DNA Viewer", img)
    key = cv2.waitKey(1) & 0xFF
    if key == ord('h'):
        show_guide = not show_guide
    elif key in [ord('q'), 27]:
        break

cap.release()
cv2.destroyAllWindows()
