import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

# Parameters
A = 1  # Amplitude
L = 4  # Length
T = 2  # Period

# Generate data
x = np.linspace(0, L, 100)

fig, ax = plt.subplots()
lines = [
    ax.plot(x, A * np.sin(n * np.pi * x / L), label=f'Harmonic {n}')[0]
    for n in range(1, 5)
]

# Set up the plot
ax.set_title('Animated Standing Waves with Harmonics')
ax.set_xlabel('x')
ax.set_ylabel('Amplitude')
ax.set_ylim(-A, A)
ax.legend()

# Animation function
def animate(t):
    for n, line in enumerate(lines, start=1):
        y = A * np.sin(n * np.pi * x / L) * np.cos(2 * np.pi * t / T)
        line.set_ydata(y)
    return lines

# Create animation
ani = animation.FuncAnimation(fig, animate, frames=np.linspace(0, T, 200), interval=20, blit=True)

# Show the plot
plt.show()