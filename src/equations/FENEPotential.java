package equations;

/*FENE (finitely-extensible nonlinear elastic) Potential
U(r) = -0.5k r₀² ln [ 1 - (r/r₀)² ]*/

import core.Monomer;
import core.Vector;

import java.util.ArrayList;

public class FENEPotential implements Potential {
    private final double rInitial;
    private final double springConst;
    private final double sigma;
    private final double epsilon;

    public FENEPotential(double sigma, double epsilon) {
        this.sigma = sigma;
        this.rInitial = 1.5 * sigma;
        this.epsilon = epsilon;
        this.springConst = (30 * epsilon) / Math.pow(sigma, 2);
    }

    @Override
    public double calcEnergy(double r) throws IllegalArgumentException {
       if (r < rInitial) {
           double rDivRInitialPow = Math.pow((r / rInitial), 2);
           double ln = Math.log(1 - rDivRInitialPow);
           double result = -0.5 * springConst * Math.pow(rInitial, 2) * ln;
           return result;
       } else {
           throw new IllegalArgumentException("r is greater than rInitial");
       }
    }

    @Override
    public String printEnergy(double result){
        return "The FENE potential energy between these two particles is: " + result + " Joules.";
    }

    @Override
    public void applyForces(ArrayList<Monomer> monomers) {
        for (int i = 0; i < monomers.size() - 1; i++) {
            Monomer m1 = monomers.get(i);
            Monomer m2 = monomers.get(i + 1);
            Vector r12 = m2.getPosition().sub(m1.getPosition());
            double r = r12.magnitude();
            if (r > 0 && r < this.rInitial) {
                double rDivR0_sq = Math.pow(r / this.rInitial, 2);
                double forceMagnitude = (-this.springConst * r) / (1.0 - rDivR0_sq);
                Vector forceDirection = r12.normalize();
                Vector forceVector = forceDirection.scale(forceMagnitude);
                m1.addForce(forceVector);
                m2.addForce(forceVector.scale(-0.1));
            }
        }
    }
}
