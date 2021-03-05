package comp559.fluid;

import java.util.LinkedList;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2f;
import javax.vecmath.Tuple2f;

import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;

/**
 * Implementation of a thin filament which is advected by fluid
 * @author kry
 */
public class Filament {
    
    private Fluid fluid;
    
    private LinkedList<Point2f> L = new LinkedList<Point2f>();

    private LinkedList<Point2f> Ltmp = new LinkedList<Point2f>();

    private double age = 0;

    private double maxAge = 60;        
    
    private static DoubleParameter defaultMaxAge = new DoubleParameter( "max age", 60, 10, 150 );
    
    private static DoubleParameter refinementThreshold = new DoubleParameter( "refinement threshold", 0.5, 0.1, 1 );
    
    private static IntParameter maxParticles = new IntParameter( "max particles", 1000, 500, 1500 );

    static public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.setBorder( new TitledBorder("Filament Controls" ));
        vfp.add( defaultMaxAge.getSliderControls(false) );
        vfp.add( refinementThreshold.getSliderControls(false));
        vfp.add( maxParticles.getSliderControls() );
        return vfp.getPanel();
    }
    
    /**
     * Gets the number of points in the filament.  This is used to help monitor the complexity of this geometry
     * @return the number of points in the filament 
     */
    public int getNumPoints() {
        return L.size();
    }
    
    public LinkedList<Point2f> getPoints() {
    	return L;
    }
    
    /**
     * creates a line at the given height
     * @param f 
     * 
     * @param h
     */
    public Filament( Fluid f, float h) {
        fluid = f;        
        maxAge = defaultMaxAge.getValue();
        for (int i = 1; i <= fluid.N + 1; i++) {
            L.add(new Point2f(i * fluid.dx, h));
        }
    }

    public Filament( Fluid f, Tuple2f x, int dir ) {
        fluid = f;        
        maxAge = defaultMaxAge.getValue();
        if( dir == 0 ) {
            for (int i = 1; i <= fluid.N + 1; i++) {
                L.add(new Point2f(i * fluid.dx, x.y));
            }
        } else if ( dir == 1 ) {
            for (int i = 1; i <= fluid.N + 1; i++) {
                L.add(new Point2f( x.x, i * fluid.dx));
            }
        }
        
    }
    
    public boolean isExpired() {
        return age > 2 * maxAge;
    }

    public void refineAndAdvect() {

        float dx = fluid.dx;
        int N = fluid.N;
        float dt = fluid.timeStepSize.getFloatValue();

        double thresh = refinementThreshold.getValue();

        // we'll copy from one list to the next
        Point2f prev = null;
        for (Point2f p : L) {
            if (prev != null) {
                if (p.distance(prev) > dx * thresh) {
                    Point2f h = new Point2f();
                    h.interpolate(prev, p, 0.5f);
                    Ltmp.add(h);
                }
            }
            Ltmp.add(p);
            prev = p;
        }
        
        LinkedList<Point2f> tmp = L;
        L = Ltmp;
        Ltmp = tmp;
        Ltmp.clear();

        // advect step...
        Point2f x1 = new Point2f();
        for (Point2f x0 : L) {
            fluid.traceParticle(x0, dt, x1);
            // clamp the particles to the boundaries
            if (x1.x < dx) x1.x = dx;
            if (x1.x > (N + 1) * dx) x1.x = (N + 1) * dx;
            if (x1.y < dx) x1.y = dx;
            if (x1.y > (N + 1) * dx) x1.y = (N + 1) * dx;
            x0.set(x1);
        }
        age += dt;

        if (age < maxAge && L.size() > maxParticles.getValue())
            maxAge = age;

    }

    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();
        gl.glBegin( GL.GL_LINE_STRIP );
        gl.glColor4d(1, 1, 1, age > maxAge ? (1.0 - (age - maxAge) / maxAge) : 1);
        for (Point2f p : L) {
            gl.glVertex2f(p.x, p.y);
        }
        gl.glEnd();
    }

}
