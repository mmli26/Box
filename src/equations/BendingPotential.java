package equations;

/*Bending Potential
U(Θᵢⱼₖ) = k [ 1 - cos(Θᵢⱼₖ) ]*/

import core.Monomer;

import core.Vector;

import java.util.ArrayList;

public class BendingPotential implements Potential {
    private final double k;
    private double epsilon;

    public BendingPotential(double k, double epsilon) throws IllegalArgumentException {
        if (k < 8 * epsilon || k > 128 * epsilon){
            throw new IllegalArgumentException("k is outside of the recommended range. Keep it inside [8ε, 128ε]");
        } else {
            this.k = k;
        }
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

    @Override
    public double calcEnergy(double theta) {
            double result = k * (1 - Math.cos(theta));
            return result;
    }

    @Override
    public String printEnergy(double result){
        return "The Bending Potential energy between these two particles is: " + result + " Joules.";
    }

    @Override
    public void applyForces(ArrayList<Monomer> monomers) {
        for (int i = 0; i < monomers.size() - 2; i++) {
            Monomer m1 = monomers.get(i);
            Monomer m2 = monomers.get(i + 1);
            Monomer m3 = monomers.get(i + 2);
            Vector r21 = m1.getPosition().sub(m2.getPosition());
            Vector r23 = m3.getPosition().sub(m2.getPosition());
            double r21_mag = r21.magnitude();
            double r23_mag = r23.magnitude();
            if (r21_mag == 0 || r23_mag == 0) continue;
            double cosTheta = r21.dot(r23) / (r21_mag * r23_mag);
            cosTheta = Math.max(-1.0, Math.min(1.0, cosTheta));
            double theta = Math.acos(cosTheta);
            double forceMagnitudeComponent = this.k * Math.sin(theta);
            if (forceMagnitudeComponent == 0) continue;
            Vector force1 = r21.cross(r23).cross(r21).normalize().scale(forceMagnitudeComponent / r21_mag);
            Vector force3 = r23.cross(r21).cross(r23).normalize().scale(forceMagnitudeComponent / r23_mag);
            Vector force2 = force1.add(force3).scale(-1.0);
            m1.addForce(force1);
            m2.addForce(force2);
            m3.addForce(force3);
        }
    }
}
