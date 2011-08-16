import time, threading
import numpy
from matplotlib.pyplot import *

x = numpy.linspace(0, 10)
y = x**2

def main():
    plot(x, x)
    draw()
    time.sleep(2)
    plot(x, y)
    draw()

thread = threading.Thread()
thread.run = main

manager = get_current_fig_manager()
manager.window.after(100, thread.start)
figure(1)
show()
