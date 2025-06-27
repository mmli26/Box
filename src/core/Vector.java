package core;

public class Vector {
    //coordinates
    private final double x;
    private final double y;
    private final double z;
    public static final Vector ZERO = new Vector(0, 0, 0);

    public Vector(double x, double y, double z){
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public double getX(){
        return x;
    }

    public double getY(){
        return y;
    }

    public double getZ(){
        return z;
    }

    public double magnitude(){
        return Math.sqrt(x*x + y*y + z*z);
    }

    public double dot(Vector v){
        return x*v.x + y*v.y + z*v.z;
    }

    public Vector cross(Vector v){
        double newX = this.y * v.z - this.z * v.y;
        double newY = this.z * v.x - this.x * v.z;
        double newZ = this.x * v.y - this.y * v.x;
        return new Vector(newX, newY, newZ);
    }

    public Vector add(Vector v){
        return new Vector(x+v.x, y+v.y, z+v.z);
    }

    public Vector sub(Vector v){
        return new Vector(x-v.x, y-v.y, z-v.z);
    }

    public Vector multiplyComponents(Vector v){
        return new Vector(x*v.x, y*v.y, z*v.z);
    }

    public Vector divideComponents(Vector v){
        return new Vector(x/v.x, y/v.y, z/v.z);
    }

    public Vector normalize() {
        double mag = this.magnitude();
        if (mag == 0.0) {
            return new Vector(0, 0, 0);
        }
        return new Vector(this.x / mag, this.y / mag, this.z / mag);
}

    public Vector scale(double scalar) {
        return new Vector(this.x * scalar, this.y * scalar, this.z * scalar);
    }
}
