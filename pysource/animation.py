
from tempfile import NamedTemporaryFile
import matplotlib.pyplot as plt

VIDEO_TAG = """<video controls>
 <source src="data:video/x-m4v;base64,{0}" type="video/mp4">
 Your browser does not support the video tag.
</video>"""

def anim_to_html(anim):
    if not hasattr(anim, '_encoded_video'):
        with NamedTemporaryFile(suffix='.mp4') as f:
            anim.save(f.name, fps=24, bitrate=3000, extra_args=['-vcodec', 'libx264'])
            plt.close()
            video = open(f.name, "rb").read()
        anim._encoded_video = video.encode("base64")
    return VIDEO_TAG.format(anim._encoded_video)

def display_animation(anim):
    from IPython.display import HTML
    import IPython.display as display

    plt.close(anim._fig)
    return HTML(anim_to_html(anim))


from matplotlib import animation
animation.Animation._repr_html_ = anim_to_html

def create_animation(settings,Propagator,plot_every=1,figsize = (5,5),transform = lambda field:abs(field)**2):
    from IPython.display import HTML
    import IPython.display as display
    from finitedifferences import Solver2D
    from plot import plot
    
    solver = Propagator(settings)
    N = settings.get(settings.symbols.nz,int)/plot_every+1
    
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(figsize=figsize)
    ax = plt.axes()
    
    field = transform(solver.get_field())
    image = plot( field , ax=ax )
    last_max = [field.max()]
    
    if len(field.data.shape) == 1:
        ax.set_ylim(bottom=0)
    
    # initialization function: plot the background of each frame
    def init():
        #image.set_data( transform(solver.get_field()).data )
        return image,

    # animation function.  This is called sequentially
    def animate(i):
        field = transform(solver.get_field()).data
            
        if len(field.shape) == 1:
            image.set_data( image.get_data()[0], field )
            curr_max = field.max()
            if last_max[0] < curr_max:
                last_max[0] = curr_max
                ax.set_ylim(0,curr_max)
                plt.draw()
        else:
            image.set_data( field )
            image.autoscale()
            
        display.display("Frame %i rendered, %i%% complete." % (i,100*(i+1)/float(N)))
        display.clear_output(wait=True)
        
        for k in range(plot_every):
            solver.step()
        
        return image,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,frames=N, interval=20, blit=True)
    
    return anim
