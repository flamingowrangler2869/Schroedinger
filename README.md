# Schrödinger
Schrödinger equation simulation in C

This version only has one potential energy function, where $V(x)=0$ if $x$ is between $0$ and $1$, and $V(x)=\infty$ otherwise.

This was sloppily made over the course of around 2 months in the winter. I never learned out to structure C projects, so it's all in one file. I also wrote it on my laptop whilst on vacation, in class, at home, and at various events, so it's even more twisted as my thoughts were completely scattered throughout development.

I plan on making a much better version in the future (probably over the summer), fitted with actual project structure, preset and custom potential functions solved via LAPACK, realistic variables (e.g. femtoseconds, nanometers, electronvolts, etc), more keybinds, momentum wavefunction, and a lot more.

| Key | Effect |
| :---: | --- |
| ↑ ↓ | change speed |
| R | reset |
| Space | pause |
| W | view wavefunction probability distribution (color: $\mathrm{arg}(\Psi)$, height: $\vert\Psi\vert^2$)|
