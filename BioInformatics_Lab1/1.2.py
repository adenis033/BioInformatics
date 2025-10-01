text = "ACGGGCATATGCGC"

alphabet = sorted(set(text))
print (alphabet)

total = len(text)

for ch in alphabet:
    count = text.count(ch)
    percent = (count / total) * 100
    print(f"{ch}: {count} ({percent}%)")