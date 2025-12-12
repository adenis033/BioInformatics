import tkinter as tk

def nw(seq1, seq2, match, mismatch, gap):
    n = len(seq1)
    m = len(seq2)

    M = [[0]*(m+1) for _ in range(n+1)]
    T = [[""]*(m+1) for _ in range(n+1)]

    for i in range(1, n+1):
        M[i][0] = i * gap
        T[i][0] = "U"

    for j in range(1, m+1):
        M[0][j] = j * gap
        T[0][j] = "L"

    for i in range(1, n+1):
        for j in range(1, m+1):
            d = M[i-1][j-1] + (match if seq1[i-1]==seq2[j-1] else mismatch)
            u = M[i-1][j] + gap
            l = M[i][j-1] + gap
            best = max(d,u,l)
            M[i][j] = best
            T[i][j] = "D" if best==d else ("U" if best==u else "L")

    a1 = ""
    a2 = ""
    mid = ""
    path = []
    i, j = n, m
    matches = 0

    while i>0 or j>0:
        path.append((i,j))
        if T[i][j] == "D":
            c1 = seq1[i-1]
            c2 = seq2[j-1]
            a1 = c1 + a1
            a2 = c2 + a2
            if c1 == c2:
                mid = "|" + mid
                matches += 1
            else:
                mid = " " + mid
            i -= 1
            j -= 1
        elif T[i][j] == "U":
            a1 = seq1[i-1] + a1
            a2 = "-" + a2
            mid = " " + mid
            i -= 1
        else:
            a1 = "-" + a1
            a2 = seq2[j-1] + a2
            mid = " " + mid
            j -= 1

    path.append((0,0))
    path.reverse()

    sim = round(matches / len(a1) * 100, 2)

    return a1, mid, a2, matches, sim, M, path


class App:
    def __init__(self, r):
        self.r = r
        r.title("Needleman–Wunsch Alignment")

        left = tk.Frame(r)
        left.grid(row=0, column=0, sticky="n")

        tk.Label(left, text="Seq1").pack()
        self.s1 = tk.Entry(left, width=30)
        self.s1.pack()
        self.s1.insert(0, "ACCGTGAAGCCAATAC")

        tk.Label(left, text="Seq2").pack()
        self.s2 = tk.Entry(left, width=30)
        self.s2.pack()
        self.s2.insert(0, "AGCGTGCAGCCAATAC")

        tk.Label(left, text="Match").pack()
        self.match = tk.Entry(left, width=5)
        self.match.insert(0, "1")
        self.match.pack()

        tk.Label(left, text="Mismatch").pack()
        self.mismatch = tk.Entry(left, width=5)
        self.mismatch.insert(0, "-1")
        self.mismatch.pack()

        tk.Label(left, text="Gap").pack()
        self.gap = tk.Entry(left, width=5)
        self.gap.insert(0, "-1")
        self.gap.pack()

        tk.Button(left, text="ALIGN", command=self.run).pack(pady=10)

        self.mat = tk.Canvas(r, width=350, height=350, bg="white")
        self.mat.grid(row=0, column=1)

        self.trace = tk.Canvas(r, width=350, height=350, bg="white")
        self.trace.grid(row=0, column=2)

        self.out = tk.Text(r, width=100, height=10)
        self.out.grid(row=1, column=0, columnspan=3, pady=10)

    def draw_matrix(self, M):
        self.mat.delete("all")
        rows = len(M)
        cols = len(M[0])
        cw = 350/cols
        ch = 350/rows
        mn = min(min(r) for r in M)
        mx = max(max(r) for r in M)
        for i in range(rows):
            for j in range(cols):
                t = (M[i][j] - mn) / (mx - mn + 1e-6)
                r = int(255 * t)
                b = int(255 * (1 - t))
                col = f"#{r:02x}00{b:02x}"
                x1 = j * cw
                y1 = i * ch
                self.mat.create_rectangle(x1, y1, x1+cw, y1+ch, fill=col, outline="")

    def draw_traceback(self, rows, cols, path):
        self.trace.delete("all")
        cw = 350/cols
        ch = 350/rows

        for i in range(rows):
            for j in range(cols):
                x1 = j*cw
                y1 = i*ch
                self.trace.create_rectangle(x1,y1,x1+cw,y1+ch,fill="#ffffdd",outline="black")

        for i,j in path:
            x1 = j*cw
            y1 = i*ch
            self.trace.create_rectangle(x1,y1,x1+cw,y1+ch,fill="#cc0000",outline="black")

    def run(self):
        s1 = self.s1.get().upper()
        s2 = self.s2.get().upper()
        match = int(self.match.get())
        mismatch = int(self.mismatch.get())
        gap = int(self.gap.get())

        a1, mid, a2, matches, sim, M, path = nw(s1, s2, match, mismatch, gap)

        self.out.delete("1.0", tk.END)
        self.out.insert(tk.END, "Show Alignment:\n")
        self.out.insert(tk.END, a1 + "\n")
        self.out.insert(tk.END, mid + "\n")
        self.out.insert(tk.END, a2 + "\n\n")
        self.out.insert(tk.END, f"Matches = {matches}\n")
        self.out.insert(tk.END, f"Length = {len(a1)}\n")
        self.out.insert(tk.END, f"Similarity = {sim}%\n")
        self.out.insert(tk.END, f"Tracing back: M[{len(M)-1},{len(M[0])-1}]\n")

        self.draw_matrix(M)
        self.draw_traceback(len(M), len(M[0]), path)


root = tk.Tk()
App(root)
root.mainloop()