from numpy import pi, sin
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import scipy.special

def signal(amp, freq):
    return amp * sin(2 * pi * freq * t)

def funcion_ajuste_FE(alpha,sigma,l):
    """Devuelve la funcion de ajuste con los parametros proporcionados"""
    a = -l
    b = l
    c = -l
    d = l

    cte = np.pi*alpha*sigma**2/(2*(b-a)*(d-c))
    raiz_2 = np.sqrt(2)
    return lambda x,y: cte*( scipy.special.erf((b-x)/raiz_2) - scipy.special.erf((a-x)/raiz_2) )*( scipy.special.erf((d-y)/raiz_2) - scipy.special.erf((c-y)/raiz_2) )


axis_color = 'lightgoldenrodyellow'
fig = plt.figure()

# Draw the plot
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.25, bottom=0.25)

limite_fantoma = 30 #cm
bins = 121
x = np.linspace(-limite_fantoma,limite_fantoma,bins)

alpha_0 = 0.00000000005
sigma_0 = 1000
l = 0.0000001
[line] = ax.plot(x,funcion_ajuste_FE(alpha_0,sigma_0,l)(x,0), linewidth=2, color='red')
ax.set_xlim([-30, 30])
ax.set_ylim([0, 0.0001])

# Add two sliders for tweaking the parameters
alpha_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axis_color)
alpha_slider = Slider(alpha_slider_ax, 'Alpha', 1E-10, 1E-9, valinit=alpha_0)
sigma_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axis_color)
sigma_slider = Slider(sigma_slider_ax, 'Sigma', 1000, 20000.0, valinit=sigma_0)

def sliders_on_changed(val):
    line.set_ydata(funcion_ajuste_FE(alpha_slider.val, sigma_slider.val, l)(x,0))
    fig.canvas.draw_idle()
alpha_slider.on_changed(sliders_on_changed)
sigma_slider.on_changed(sliders_on_changed)

# Add a button for resetting the parameters
reset_button_ax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
reset_button = Button(reset_button_ax, 'Reset', color=axis_color, hovercolor='0.975')
def reset_button_on_clicked(mouse_event):
    sigma_slider.reset()
    alpha_slider.reset()
reset_button.on_clicked(reset_button_on_clicked)

# Add a set of radio buttons for changing color
color_radios_ax = fig.add_axes([0.025, 0.5, 0.15, 0.15], axisbg=axis_color)
color_radios = RadioButtons(color_radios_ax, ('red', 'blue', 'green'), active=0)
def color_radios_on_clicked(label):
    line.set_color(label)
    fig.canvas.draw_idle()
color_radios.on_clicked(color_radios_on_clicked)

plt.show()