package equations;

import core.Monomer;
import core.Vector;
// NEEDS WORK WITH <. . .>
/*Persistence Length
 lp = -lb / ln(<cos(θᵢⱼₖ)>)*/
public class PersistenceLength {
    private final double lb;

    public PersistenceLength(double lb){
        this.lb = lb;
    }

    public double setTheta(Monomer m1, Monomer m2, Monomer m3){
        Vector pos1 = m1.getPosition();
        Vector pos2 = m2.getPosition();
        Vector pos3 = m3.getPosition();
        Vector r21 = pos1.sub(pos2);
        Vector r23 = pos3.sub(pos2);
        double dotProduct = r21.dot(r23);
        double magProduct = r21.magnitude() * r23.magnitude();
        if (magProduct == 0.0) {
            return Math.PI;
        }
        double cosTheta = Math.max(-1.0, Math.min(1.0, dotProduct / magProduct));
        double theta = Math.acos(cosTheta);
        return theta;
    }

    public double calcPL(double lb, double theta){
        double denom = Math.log(Math.cos(theta));
        double result = -lb / denom;
        return result;
    }

    public String printPL(double result){
        return "The Persistence Length is: " + result + "Å";
    }
}
