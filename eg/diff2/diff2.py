with open("diff2.dat", "w") as f:
    for x1 in range(10):
        for x2 in range(10):
            if x1 < x2: break
            f.write(f"{x1 - x2} {x1} {x2}\n")
