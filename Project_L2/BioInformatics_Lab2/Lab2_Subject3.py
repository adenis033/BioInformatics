# Design an app using the A.I., with contains a GUI that allows the user
# to choose a FASTA file

# The content of the file should be analyzed by using a sliding window
# of 30 positions. The content for each sliding window
# should be used in order to extract the relative freq
# of the symbols found in the alphabet of the sequence, thus your input
# should be the DNA seq from the FASTA file and the output should be the
# values of the relative freq of each symbol from the sequence. Translate
# in lines on a chart. Thus your chart in the case of DNA should have
# 4 lines which reflect the values found over the sequence

# Design an application using AI, which contains a GUI that allows the user to select a FASTA file.
# The content of the file should be analyzed by using a sliding window of 30 positions.
# The content of each sliding window should be used in order to extract the relative frequencies of the symbols found in the alphabet of the sequence
# Thus, your input should be the DNA sequence from the FASTA file, and the output should be values of the relative frequencies for each symbol from the alphabet (sequence).
# Translate in lines on a chart, thus your chart in the case of DNA should contain 4 lines, one for each symbol found over the sequence.

import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from Bio import SeqIO
import numpy as np

DNA_ALPHABET = ['A', 'C', 'G', 'T']
WINDOW_SIZE = 30

class FastaAnalyzerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FASTA Nucleotide Frequency Analyzer")
        self.root.geometry("1000x750")

        self.current_file_path = None

        top_frame = tk.Frame(self.root, pady=5)
        top_frame.pack(side="top", fill="x")

        chart_frame = tk.Frame(self.root)
        chart_frame.pack(side="top", fill="both", expand=True, padx=10, pady=10)

        # Button to select FASTA file
        self.btn_select_file = tk.Button(
            top_frame,
            text="Select FASTA File",
            command=self.select_file,
            font=("Helvetica", 12)
        )
        self.btn_select_file.pack(side="left", padx=20)

        # Label to show selected file path
        self.lbl_file_path = tk.Label(
            top_frame,
            text="No file selected",
            font=("Helvetica", 10),
            fg="gray"
        )
        self.lbl_file_path.pack(side="left")

        # Setup matplotlib Figure and Canvas
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=chart_frame)
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)

        self.setup_initial_chart()

    def setup_initial_chart(self):
        self.ax.clear()
        self.ax.set_title("Sliding Window Nucleotide Frequencies")
        self.ax.set_xlabel("Window Start Position")
        self.ax.set_ylabel("Relative Frequency (%)")  # Y-axis in percent
        self.ax.grid(True, linestyle='--', alpha=0.6)
        self.ax.text(0.5, 0.5, "Please select a FASTA file to begin analysis.",
                     horizontalalignment='center', verticalalignment='center',
                     transform=self.ax.transAxes, fontsize=14, color='gray')
        self.canvas.draw()

    def select_file(self):
        file_path = filedialog.askopenfilename(
            title="Select a FASTA file",
            filetypes=(("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*"))
        )
        if file_path:
            self.current_file_path = file_path
            self.lbl_file_path.config(text=file_path, fg="black")
            self.process_and_plot()

    def apply_smoothing(self, data, window_size):
        if window_size <= 1:
            return np.array(data)
        kernel = np.ones(window_size) / window_size
        smoothed_data = np.convolve(data, kernel, mode='valid')
        return smoothed_data

    def process_and_plot(self):
        if not self.current_file_path:
            return

        try:
            fasta_sequences = SeqIO.parse(open(self.current_file_path), 'fasta')
            first_record = next(fasta_sequences)
            sequence = str(first_record.seq).upper()

            if len(sequence) < WINDOW_SIZE:
                messagebox.showerror("Error", "Sequence length is too short.")
                return

            frequencies, window_positions = self.sliding_window_analysis(sequence)

            # Always smoothing
            smoothing_window = 51  # Fixed smoothing window
            smoothed_frequencies = {}
            for symbol in DNA_ALPHABET:
                # Apply smoothing
                smoothed_frequencies[symbol] = self.apply_smoothing(
                    np.array(frequencies[symbol]) * 100,  # Convert to percentages
                    smoothing_window
                )

            offset = (smoothing_window - 1) // 2
            adjusted_positions = window_positions[offset:-offset]

            self.plot_frequencies(
                smoothed_frequencies,
                adjusted_positions,
                first_record.id,
                smoothing_window
            )

        except StopIteration:
            messagebox.showerror("Error", "The selected file is empty or not a valid FASTA file.")
        except Exception as e:
            messagebox.showerror("An Error Occurred", f"An unexpected error occurred: {e}")
            self.setup_initial_chart()

    def sliding_window_analysis(self, sequence):
        frequencies = {symbol: [] for symbol in DNA_ALPHABET}
        window_positions = []
        for i in range(len(sequence) - WINDOW_SIZE + 1):
            window = sequence[i : i + WINDOW_SIZE]
            for symbol in DNA_ALPHABET:
                count = window.count(symbol)
                relative_freq = count / WINDOW_SIZE
                frequencies[symbol].append(relative_freq)
            window_positions.append(i)
        return frequencies, window_positions

    def plot_frequencies(self, frequencies, positions, sequence_id, smoothing_window=0):
        self.ax.clear()
        colors = {'A': 'blue', 'C': 'green', 'G': 'orange', 'T': 'red'}

        for symbol in DNA_ALPHABET:
            self.ax.plot(positions, frequencies[symbol], label=f"Freq of '{symbol}'", color=colors[symbol])

        title = f"Nucleotide Frequencies for: {sequence_id} (Smoothed)"
        self.ax.set_title(title)
        self.ax.set_xlabel("Window Start Position")
        self.ax.set_ylabel("Relative Frequency (%)")
        self.ax.set_ylim(0, 100)  # Percent range
        self.ax.legend()
        self.ax.grid(True, linestyle='--', alpha=0.6)
        self.canvas.draw()

def main():
    root = tk.Tk()
    app = FastaAnalyzerApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()