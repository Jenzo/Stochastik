package test;

import java.math.RoundingMode;
import java.text.DecimalFormat;

import util.Punkt;
import util.Stochastik;
import util.Stochastik.BernoulliTuple;

public class Start {

	final static DecimalFormat df = new DecimalFormat("#.#####");

	public static void main(String[] args) {
		df.setRoundingMode(RoundingMode.HALF_UP);

		double p;
		int n = 60;
		int k = 5;
		System.out.println(n + " über " + k + ": " + Stochastik.nUeberK(n, k));

		System.out.println(Stochastik.BernulliTuple(100, 0.51, 50));
		System.out.println(Stochastik.BernulliAbschätzung(100, 0.51, 50));

		// Aufgabe 35 c.i) Brauerei
		// 500(n) Kisten werden ausgeliefert, Wahrscheinlichkeit p, dass 1 Kiste zu
		// Bruch geht(p) 0.002
		// Wahrscheinlichkeit, dass 3 Kisten (k) zu Bruch gehen ?
		// n ist groß ; p ist klein ; k nicht so groß
		// ==> Poisson-Nährung für Binominalverteilung
		n = 500;
		p = 0.002;
		k = 3;

		BernoulliTuple b = new BernoulliTuple(n, p, k);
		double result = Stochastik.PoissonNährung(b);

		System.out.println("Aufgabe 35 c.i) -> " + result);

		// c.ii) weniger als 3 Kisten zu Bruch gehen
		// Ki = {0,1,2}

		result = Stochastik.PoissonNährung(n, p, 3);

		System.out.println("Aufgabe 35 c.ii) -> " + result);

		// Aufgabe 36
		// 100 Bernoulli-Experimente werden durchgeführt mit einer Wahrscheinlichkeit
		// von von p = 0.8 für ein bestimmtes Ereignis.
		// i) Wahrscheinlichkeit, dass Ereignis nicht weniger als 75 und nicht mehr als
		// 90 mal eintritt
		// k1 = 74 , k = 75 , k2 = 90
		// n*p*q > 9 ? -> n = 100, p = 0.8, q = 1 - 0.8 = 0.2
		// n*p*q = 16
		// sigma = sqrt(n*p*q)
		// somit T(k1 < k <= k2) = T((k2 - µ + 0.5)/sigma) - T((k1 - µ + 0.5)/sigma))

		result = Stochastik.x(100, 0.8, 90, 74);
		System.out.println("Aufgabe 36 -> " + result);

		System.out.println(Stochastik.bla(0));

		double x = 0.5;
		System.out.println("Gaußsche Normalverteilung für x = " + x + " : "
				+ df.format(Stochastik.GaußscheNormalverteilungTabellenwert(x)));

		// Aufgabe 58) Regression mit linearer Funktion
		Punkt[] punkte = new Punkt[] { new Punkt(1, 3.4), new Punkt(1.5, 5.2), new Punkt(2, 7), new Punkt(2.5, 8.3),
				new Punkt(3, 10.5), new Punkt(3.5, 12), new Punkt(4, 13.7) };

		System.out.println("Aufgabe 58 -> " + Stochastik.getProportionaleFunktion(punkte));

		punkte = new Punkt[] { new Punkt(2, 3), new Punkt(5, 4), new Punkt(8, 8), new Punkt(10, 11) };
		System.out.println("Aufgabe 59 -> " + Stochastik.getLineareFuntion(punkte));

		punkte = new Punkt[]{new Punkt(0, 1.6), new Punkt(1, 1.8), new Punkt(1.5, 2), new Punkt(2, 2.2), new Punkt(2.5, 2.5)};
		System.out.println("Aufgabe 60 -> " + Stochastik.getExponentialFunktion(punkte));

		punkte = new Punkt[]{new Punkt(1, 0.6), new Punkt(1.5, 0.9), new Punkt(2, 1.4), new Punkt(2.5, 2), new Punkt(3, 2.6)};
		System.out.println("Aufgabe 61 -> " + Stochastik.getPotenzFunktion(punkte));

	}

}
