package core;

public class Monomer {
    private Vector position;
    private Vector velocity;
    private Vector force;
    private final double mass;
    private final int type;
    private final int idx;

    public Monomer(Vector position, Vector velocity, Vector force, double mass, int type, int idx){
        this.position = position;
        this.velocity = velocity;
        this.force = force;
        this.mass = mass;
        this.type = type;
        this.idx = idx;
    }

    public Vector getPosition(){
        return position;
    }

    public Vector getVelocity(){
        return velocity;
    }

    public Vector getForce(){
        return force;
    }

    public double getMass(){
        return mass;
    }

    public int getType(){
        return type;
    }

    public int getIdx(){
        return idx;
    }

    public void setPosition(Vector position){
        this.position = position;
    }

    public void setVelocity(Vector velocity){
        this.velocity = velocity;
    }

    private void setForce(Vector force){
        this.force = force;
    }

    public void clearForce(){
        this.force = Vector.ZERO;
    }

    public void addForce(Vector f){
        this.force = this.force.add(f);
    }
}
