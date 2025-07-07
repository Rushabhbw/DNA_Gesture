# DNA_Gesture
Hand-Gesture-Controlled DNA Sequence Viewer and Analyzer
This project allows interactive exploration and analysis of DNA sequences using hand gestures via webcam. It integrates OpenCV, MediaPipe, and BioPython to visualize DNA chunks, translate them into amino acids, identify ORFs, and generate GC content and codon frequency plots—all controllable with finger gestures.

Folder Structure:
dna_gesture/
├── codon_bar.png            # Bar plot of top 10 codons (auto-generated)
├── dna_gesture_viewer.py    # Main Python script with gesture controls
├── gc_plot.png              # GC% content plot across the sequence (auto-generated)
├── gene_report.html         # HTML report with DNA chunk, AA, ORFs, plots (auto-generated)
├── gene_summary.txt         # Text summary with ORFs and basic stats (auto-generated)
└── your_sequence.fasta      # Example FASTA file to test the system

🚀 Features
✅ Hand Gesture Navigation:
🖐️ 1 finger → Next DNA chunk
🖐️ 2 fingers → Previous DNA chunk
🖐️ 3 fingers → Show GC% of current chunk
🖐️ 4 fingers → Export HTML report (gene_report.html)
🖐️ 5 fingers → Save text summary (gene_summary.txt)
> Press 'h' to show/hide gesture guide
> Press 'q' or ESC to exit.
✅ Visual DNA Chunk Viewer:
Color-coded base display:
A: Blue, T: Pink, G: Green, C: Purple.
> Amino acid translation of current chunk displayed.
> Top detected ORFs with positions and translations shown.
> GC% and codon frequency plots displayed live.
✅ Automated Analysis and Reports:
> gc_plot.png: GC% plot across the entire sequence.
> codon_bar.png: Top 10 codon usage chart.
> gene_report.html: Clean, detailed HTML report.
> gene_summary.txt: Quick report for notes or documentation.

⚙️ Requirements
Python 3.8+, Webcam
Dependencies:
opencv-python, mediapipe, biopython, matplotlib, numpy, tkinter (preinstalled with Python).
Install all with:
pip install opencv-python mediapipe biopython matplotlib numpy

🧪 Usage
1️⃣ Place your .fasta DNA sequence file in the folder.
2️⃣ Run:
> python dna_gesture_viewer.py
3️⃣ A file dialog will appear to select your .fasta file.
4️⃣ Use hand gestures in front of your webcam to:
Navigate chunks, View GC%, Export HTML reports, Generate text summaries, and explore your DNA sequence hands-free.

📊 Outputs
> gc_plot.png – GC% plot across your sequence.
> codon_bar.png – Top 10 codon frequency chart.
> gene_report.html – Shareable, clean HTML report.
> gene_summary.txt – Quick summary with ORFs, GC%, and translation.

💡 Applications
✅ Teaching bioinformatics interactively
✅ Hands-free DNA analysis while working at the bench
✅ Gesture-controlled UI experiments for scientific data
✅ Quick DNA chunk analysis and ORF finding

✍️ Author
Rushabh Wakade

📜 License
This project is open for academic, learning, and personal use. For any publication or commercial use, please credit the author.
