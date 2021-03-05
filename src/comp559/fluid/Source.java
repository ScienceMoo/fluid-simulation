package comp559.fluid;

import javax.vecmath.Point2f;

/**
 * Heating or cooling source in the fluid
 * 
 * @author kry
 */
public class Source {
    
    /** position of the heat source */
    Point2f location = new Point2f();
    
    /** heating or cooling rate */
    public float amount = 0;
    
    /** flag to denote that the mouse is over the source */
    boolean highlight = false;
    
    /**
     * Creates a source with the given amount
     * @param p
     * @param a
     */
    public Source( Point2f p, float a ) {
        location.set(p);
        amount = a;
        highlight = false;
    }
    
}