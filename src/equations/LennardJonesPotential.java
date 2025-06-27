package equations;

/*Lennard-Jones Potential
U(r) = 4ε [ (σ/r)¹² - (σ/r)⁶ ]*/

import core.Monomer;
import core.Vector;

import java.util.ArrayList;

public class LennardJonesPotential implements Potential {
    private final double epsilon;
    private final double sigma;
    private final double rCutoff;

    public LennardJonesPotential(double epsilon, double sigma) {
        this.epsilon = epsilon;
        this.sigma = sigma;
        this.rCutoff = 2.5 * sigma;
    }

    @Override
    public double calcEnergy(double r){
        if (r <= rCutoff) {
            double sigmaOverR = sigma / r;
            double repulsiveTerm = Math.pow(sigmaOverR, 12);
            double attractiveTerm = Math.pow(sigmaOverR, 6);
            double result = 4.0 * epsilon * (repulsiveTerm - attractiveTerm);
            return result;
        } else {
            return 0.0;
        }
    }

    @Override
    public String printEnergy(double result){
        return "The Lennard-Jones potential energy between these two particles is: " + result + " Joules.";
    }

    @Override
    public void applyForces(ArrayList<Monomer> monomers) {
        for (int i = 0; i < monomers.size(); i++) {
            for (int j = i + 1; j < monomers.size(); j++) {
                if (Math.abs(i - j) > 1) {
                    Monomer m1 = monomers.get(i);
                    Monomer m2 = monomers.get(j);
                    Vector r12 = m2.getPosition().sub(m1.getPosition());
                    double r = r12.magnitude();
                    if (r > 0) {
                        double sigmaOverR = this.sigma / r;
                        double attractive = Math.pow(sigmaOverR, 6);
                        double repulsive = Math.pow(sigmaOverR, 12);
                        double forceMagnitude = (24.0 * this.epsilon / r) * (2.0 * repulsive - attractive);
                        Vector forceDirection = r12.normalize();
                        Vector forceVector = forceDirection.scale(forceMagnitude);
                        m1.addForce(forceVector);
                        m2.addForce(forceVector.scale(-1.0));
                    }
                }
            }
        }
    }
}
