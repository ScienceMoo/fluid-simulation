package comp559.fluid;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2f;
import javax.vecmath.Tuple2f;
import javax.vecmath.Vector2f;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import com.jogamp.opengl.util.gl2.GLUT;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.EasyViewer;
import mintools.viewer.Interactor;
import mintools.viewer.SceneGraphNode;


/**
 * Fluid app, provides drawing code, keyboard and mouse controls for interaction 
 * and provides a UI panel.
 * @author kry
 */     
public class FluidApp implements SceneGraphNode, Interactor {

    /** 
     * starts the application 
     * @param args 
     */
    public static void main( String[] args ) {
        new FluidApp();
    }
            
    private Fluid fluid = new Fluid();
    
    private EasyViewer ev;
    
    /**
     * Creates a fluid simulation assignment application
     */
    public FluidApp() {
        fluid.setup();
        ev = new EasyViewer( "Stable Fluid", this, new Dimension( 512, 512 ), new Dimension(600,600)  );
        ev.addInteractor(this);
    }
   
    private List<Filament> loops = new LinkedList<Filament>();
    
    private Point2f Xdrag = null;
    
    private Point2f X = new Point2f();
        
    private Point2f X1 = new Point2f();
    
    private Point2f X0 = new Point2f();
        
    boolean b1, b3;
    
    @Override
    public void attach(Component component) {
        component.addMouseListener( new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                if ( e.getButton() == MouseEvent.BUTTON1 ) {
                    for ( Source s : fluid.sources ) {
                        if (s.highlight) return;
                    }
                    setPosFromMouse(e, X);                    
                    Source s = new Source( X, 0 );
                    s.highlight = true;
                    fluid.sources.add( s );
                }
                if ( e.getButton() == MouseEvent.BUTTON2 ) {
                    setPosFromMouse(e, X);
                    loops.add( new Filament( fluid, X, 0 ) );
                    loops.add( new Filament( fluid, X, 1 ) );
                }
                if ( e.getButton() == MouseEvent.BUTTON3 ) {
                    Source toremove = null;
                    for ( Source s : fluid.sources ) {                        
                        if (s.highlight) {
                            toremove = s;
                            break;
                        }
                    }
                    if ( toremove != null ) fluid.sources.remove(toremove);
                }
            }
            @Override
            public void mousePressed(MouseEvent e) {
                if ( e.getButton() == MouseEvent.BUTTON1 ) b1 = true;                                
                if ( e.getButton() == MouseEvent.BUTTON3 ) {
                    setPosFromMouse(e, X1 );
                    X0.set( X1 );
                    b3 = true;
                }
            }
            @Override
            public void mouseReleased(MouseEvent e) {
                if ( e.getButton() == MouseEvent.BUTTON1 ) {
                    b1 = false;
                    Xdrag = null;
                }
                if ( e.getButton() == MouseEvent.BUTTON3 ) {
                    b3 = false;                
                    setPosFromMouse(e, X1 );
                    X0.set( X1 );
                }
            }
        });
        component.addMouseWheelListener( new MouseWheelListener() {
            @Override
            public void mouseWheelMoved(MouseWheelEvent e) {                 
                 int ticks = -e.getWheelRotation();
                 for ( Source s : fluid.sources ) {
                     if ( s.highlight ) {
                         s.amount += ticks / 10.0;
                     }
                 }
            }
        });
        component.addMouseMotionListener( new MouseMotionListener() {
            @Override
            public void mouseDragged(MouseEvent e) {
                
                if ( b1  && Xdrag != null ) {
                    setPosFromMouse( e, Xdrag );
                }
                if ( b3 ) {                    
                    setPosFromMouse( e, X1 );                    
                }                
            }
            @Override
            public void mouseMoved(MouseEvent e) {
                setPosFromMouse(e, X);                
                float min = Float.MAX_VALUE;
                Source closest = null;
                for ( Source s : fluid.sources ) {
                    s.highlight = false;
                    float d = s.location.distance(X);
                    if ( d < 5/scale  && d < min ) {
                        min = d;
                        closest = s;                            
                    } 
                }                    
                if ( closest != null ) {
                    Xdrag = closest.location;
                    closest.highlight = true;
                }                
            }
        });
        component.addKeyListener( new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent e) {
                if ( e.getKeyCode() == KeyEvent.VK_S ) {
                    stepRequest = true;
                } else if ( e.getKeyCode() == KeyEvent.VK_R ) {
                    resetRequest = true;
                } else if ( e.getKeyCode() == KeyEvent.VK_C ) {
                    createFilamentsRequest = true;
                } else  if ( e.getKeyCode() == KeyEvent.VK_SPACE ) {
                    run.setValue( !run.getValue() );
                } else if ( e.getKeyCode() == KeyEvent.VK_ENTER ) {
                    // toggle recording of steps to png files
                    record.setValue( ! record.getValue() );
                } else if ( e.getKeyCode() == KeyEvent.VK_UP ) {
                	for ( Source s : fluid.sources ) {
                        if ( s.highlight ) {
                            s.amount += 1 / 10.0;
                        }
                    }
                } else if ( e.getKeyCode() == KeyEvent.VK_DOWN ) {
                	for ( Source s : fluid.sources ) {
                        if ( s.highlight ) {
                            s.amount -= 1 / 10.0;
                        }
                    }
                }
            }            
        });
    }
    
    /**
     * Computes a position in fluid coordinates from the mouse position
     * @param e
     * @param x
     */
    private void setPosFromMouse( MouseEvent e, Tuple2f x ) {
        //TODO: fix DPI stuff
        int mx = 2 * e.getX();
        int my = 2 * e.getY();
        x.x = (float) (( mx - offset ) / scale);
        x.y = (float) (( my - offset ) / scale);
    }
    
    private double scale = 30;
    
    /** border size for empty space around the fluid grid display */
    private double offset = 30;
    
    private boolean stepRequest = false;
    
    private boolean resetRequest = false;
           
    private boolean createFilamentsRequest = false;
                
    @Override
    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();
        
        // resetting the viewing scale is only necessary on a resize, but we'll do it every time.
        setViewingScale(drawable);
        
        // set mouse forces in the fluid for interaction
        fluid.setMouseMotionPos( X0, X1 );
        X0.set( X1 );
        
        if ( resetRequest ) {
            loops.clear();
            fluid.setup();
            resetRequest = false;
        }        
        if ( createFilamentsRequest ) {
            for ( Source s : fluid.sources ) {
                loops.add( new Filament( fluid, s.location.y ) );
            }
            createFilamentsRequest = false;
        }
        if ( stepRequest || run.getValue() ) {
            fluid.step();
            LinkedList<Filament> toremove = new LinkedList<Filament>();
            for ( Filament l : loops ) {
                l.refineAndAdvect();
                if ( l.isExpired() ) {
                    toremove.add( l );
                }
            }        
            loops.removeAll( toremove ); 
        }        
        
        EasyViewer.beginOverlay(drawable);
        
        gl.glPushMatrix();
        gl.glTranslated( offset, offset, 0 );
        gl.glScaled( scale, scale, scale );
                
    	int N = fluid.N;
        
        // draw the scalar field
        float[] S = fluid.temperature0;        
        float dx = fluid.dx;
        Vector2f x = new Vector2f();        
        float cs = colorScale.getFloatValue();
        if ( drawScalars.getValue() ) {
        	int rf = refineFactor.getValue(); //  refine the grid by a factor of 4 ?
            int R = (N+2) * rf;
            int low = drawBoundaryCells.getValue() ? 0 : rf;
            int high = (drawBoundaryCells.getValue() ? N+2 : N+1) * rf;
            if ( drawSmooth.getValue() ) {
            	if ( highQualitySmooth.getValue() ) {
	                for ( int i = low; i < high; i++ ) {                    
	                    gl.glBegin( GL2.GL_QUAD_STRIP );                    
	                    for ( int j = low; j <= high; j++ ) {
	                        x.x = ((float)i)/R * (N+2) * dx;
	                        x.y = ((float)j)/R * (N+2) * dx;
	                        
	                        float s = fluid.interpolate(x, S) * cs;
	                       
	                        gl.glColor3d( s>0?s:0, 0.125, s<0?-s:0 );
	                        gl.glVertex2d( x.x, x.y );
	                        x.x = ((float)(i+1))/R * (N+2) * dx;
	                        s = fluid.interpolate(x, S) * cs;
	                        gl.glColor3d( s>0?s:0, 0.125 ,s<0?-s:0 );
	                        gl.glVertex2d( x.x, x.y );
	                    }
	                    gl.glEnd();
	                }
            	} else {
	                for ( int i = 0; i <= N; i++ ) {      
	                    gl.glBegin( GL2.GL_QUAD_STRIP );                    
	                    for ( int j = 0; j <= N+1; j++ ) {
	                        double s = S[fluid.IX(i,j)] * cs;
	                        gl.glColor3d( s>0?s:0, 0.125, s<0?-s:0 );
	                        x.x = (i + 0.5f) * dx;
	                        x.y = (j + 0.5f) * dx;
	                        gl.glVertex2d( x.x, x.y );           
	                        
	                        s = S[fluid.IX(i+1,j)] * cs;
	                        gl.glColor3d( s>0?s:0, 0.125, s<0?-s:0 );
	                        x.x = (i+1 + 0.5f) * dx;
	                        x.y = (j + 0.5f) * dx;
	                        gl.glVertex2d( x.x, x.y );                        
	                    }  
	                    gl.glEnd();
	                }
                }
                    
            } else {
                gl.glPointSize(3);
                gl.glBegin( GL.GL_POINTS  );
                for ( int i = 0; i <= N+1; i++ ) {                                    
                    for ( int j = 0; j <= N+1; j++ ) {
                        double s = S[fluid.IX(i,j)] * cs;
                        gl.glColor3d( s>0?s:0, 0.25, s<0?-s:0 );
                        x.x = (i + 0.5f) * dx;
                        x.y = (j + 0.5f) * dx;
                        gl.glVertex2d( x.x, x.y );                        
                    }                    
                }
                gl.glEnd();
            }
        }
        
        // draw the grid
        if ( drawGrid.getValue() ) {
            gl.glDisable( GL2.GL_LIGHTING );        
            gl.glBegin( GL.GL_LINES );
            gl.glColor3f( 0.2f, 0.2f, 0.2f );
            for ( int i = 0; i <= N+2; i++ ) {
                gl.glVertex2d( 0, i*dx );
                gl.glVertex2d( (N+2)*dx, i*dx );
                gl.glVertex2d( i*dx, 0 );
                gl.glVertex2d( i*dx, (N+2)*dx );                
            }
            gl.glEnd();
        }
        
        // draw the fluid boundary box 
        if ( drawBox.getValue() ) {
            gl.glBegin( GL.GL_LINES );
            gl.glColor3f( 1, 1, 1 );
            gl.glVertex2d( dx, 1*dx );
            gl.glVertex2d( (N+1)*dx, 1*dx );
            gl.glVertex2d( 1*dx, dx );
            gl.glVertex2d( 1*dx, (N+1)*dx );
            gl.glVertex2d( dx, (N+1)*dx );
            gl.glVertex2d( (N+1)*dx, (N+1)*dx );
            gl.glVertex2d( (N+1)*dx, dx );
            gl.glVertex2d( (N+1)*dx, (N+1)*dx );
            gl.glEnd();
        }
        
        // draw the velocities
        double vds = velocityDisplayScale.getValue();        
        if ( drawVelocities.getValue() ) {
            Vector2f pp = new Vector2f();
            Vector2f pv = new Vector2f();
            for ( int i = 0; i <= N+1; i++ ) {
                for ( int j = 0; j <= N+1; j++ ) {
                    pp.x = (i + 0.5f) * dx;
                    pp.y = (j + 0.5f) * dx;                    
                    fluid.getVelocity(pp, pv);
                    gl.glBegin( GL.GL_LINES );
                    gl.glColor4f( 0,1,0,0.5f );
                    gl.glVertex2d( pp.x, pp.y );
                    gl.glVertex2d( pp.x + pv.x * vds, pp.y + pv.y * vds );
                    gl.glEnd();
                }
            }
        }

        // draw the sources
        if ( drawSources.getValue() ) {
            final Vector2f velocity = new Vector2f();

            for ( Source s : fluid.sources ) {
                x.set( s.location );
                fluid.getVelocity(x, velocity);
        
                gl.glPointSize( s.highlight ? 10 : 5 );
                gl.glBegin( GL.GL_POINTS );
                gl.glColor4f( 0,1,0,1 );
                gl.glVertex2d( x.x, x.y );
                gl.glEnd();
                gl.glBegin( GL.GL_LINES );
                gl.glColor4f( 0,1,0,1 );
                gl.glVertex2d( x.x, x.y );
                gl.glVertex2d( x.x + velocity.x * vds, x.y + velocity.y * vds );
                gl.glEnd();
                if ( s.highlight ) {
                    gl.glColor3f(1,1,1);
                    gl.glRasterPos2d(x.x + dx*0.5, x.y - dx*0.5);        
                    EasyViewer.glut.glutBitmapString(GLUT.BITMAP_8_BY_13, "" + s.amount );
                }
            }
        }
        
        // draw the filaments
        for ( Filament l : loops ) {
            l.display(drawable);
        }
        
        gl.glPopMatrix();
        
        if ( drawInfo.getValue() ) {
            gl.glColor3f(1,1,1);
            gl.glRasterPos2d(15,15);        
            String text = "time = " + fluid.elapsed + 
            		"\nT = " + fluid.getReferenceTemperature();
            EasyViewer.printTextLines( drawable, text );
        }
        
        EasyViewer.endOverlay(drawable);
        
        if ( stepRequest || run.getValue() ) {
            if ( record.getValue() ) {
            	if ( !recordingStarted ) {
            		ev.startRecord( videoFileName.getText(), 30 ); // 30 FPS video
            		recordingStarted = true;
            		numRecordedFrames = 0;
            	}
                ev.record(drawable);
            	numRecordedFrames++;
                
                EasyViewer.beginOverlay( drawable );
                String text =  "RECORDED: " + numRecordedFrames + " frames to " + videoFileName.getText();
                gl.glDisable( GL2.GL_LIGHTING );
                gl.glColor4f( 1, 0, 0, 1 );           
                EasyViewer.printTextLines( drawable, text, 10, drawable.getSurfaceHeight()-20, 10, GLUT.BITMAP_HELVETICA_10 );
                gl.glEnable( GL2.GL_LIGHTING );
                EasyViewer.endOverlay(drawable);
            }
            stepRequest = false;
        }
        if ( !record.getValue() && recordingStarted ) {
    		ev.finishRecord();
    		recordingStarted = false;
        }
    }
    
    private BooleanParameter run = new BooleanParameter( "animate", false );
    private DoubleParameter velocityDisplayScale = new DoubleParameter( "velocity display scale", 0.1, 0.01, 100 );
    private BooleanParameter drawVelocities = new BooleanParameter( "draw velocities", true );
    private BooleanParameter drawBoundaryCells = new BooleanParameter( "draw boundary cells", true );
    private BooleanParameter drawScalars = new BooleanParameter( "draw scalar field", true );
    private BooleanParameter drawSmooth = new BooleanParameter( "draw smooth", true );
    private IntParameter refineFactor = new IntParameter( "smooth refinement factor", 1, 1, 8 );
    private BooleanParameter highQualitySmooth = new BooleanParameter( "draw smooth in high quality", true );
    private BooleanParameter drawGrid = new BooleanParameter( "draw grid", true );
    private BooleanParameter drawBox = new BooleanParameter( "draw box", true );
    private BooleanParameter drawSources = new BooleanParameter( "draw sources", true );
    private DoubleParameter colorScale = new DoubleParameter( "color scale", 1, 1e-1, 1e5 );
    private BooleanParameter drawInfo = new BooleanParameter( "draw info", true );
    private BooleanParameter record = new BooleanParameter( "record (press ENTER in canvas to toggle)", false );
    private JTextField videoFileName = new JTextField("demo.mp4");
    private boolean recordingStarted = false;    
    private int numRecordedFrames = 0;

    @Override
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.add( run.getControls() );
        vfp.add( videoFileName );
        vfp.add( record.getControls() );
        
        vfp.add( fluid.getControls() );
        
        vfp.add( Filament.getControls() );
        
        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder( new TitledBorder("display controls" ));
        vfp2.add( velocityDisplayScale.getSliderControls( true ) );
        vfp2.add( colorScale.getSliderControls(true ) );
        vfp2.add( drawVelocities.getControls() );
        vfp2.add( drawScalars.getControls() );
        vfp2.add( drawSmooth.getControls() );
        vfp2.add( refineFactor.getSliderControls() );
        vfp2.add( highQualitySmooth.getControls() );
        vfp2.add( drawGrid.getControls() );
        vfp2.add( drawBoundaryCells.getControls() );
        vfp2.add( drawBox.getControls() );
        vfp2.add( drawSources.getControls() );
        vfp2.add( drawInfo.getControls() );
        
        vfp.add( vfp2.getPanel() );
                
        return vfp.getPanel();
    }
    
    /**
     * Adjusts the scale of the display to match the window
     * @param drawable
     */
    private void setViewingScale( GLAutoDrawable drawable ) {
        width = drawable.getSurfaceWidth();
        height = drawable.getSurfaceHeight();
        double v = height;
        if ( width < height ) {
            v = width;
        }
        scale = (v - offset*2) / (fluid.dx * (fluid.N+2 ));
    }
    
    private int width;
    private int height;
    
    @Override
    public void init(GLAutoDrawable drawable) {
        setViewingScale(drawable);
    }
    
}
