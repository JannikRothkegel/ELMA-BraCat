#ifndef Histogram1D_H
#define Histogram1D_H

class Histogram1D
{
public:

	Histogram1D(double min=0.0, double max=1.0, int nrbins=10, double initValue=0.0){

		NrBins = nrbins;


		N= 0;
		Min = min;


		Max = max;


		wert = 0;


		dI = (max-min)/(1.0*nrbins);



		MaximumValue = 0;
		//Werte = new double[NrBins1+1][NrBins2+1];

		Werte = new double[NrBins+1];
		for(int i = 0; i < (NrBins+1); ++i) {
			Werte[i] = initValue;
		}
	}

	virtual ~Histogram1D(){

		delete [] Werte;

		}

	void AddValue(double x1) {

			wert = x1;
			N++;


			if( ( x1 > Max))
			{
				std::ostringstream temp;  //temp as in temporary
				temp<<x1 << " out of "<< +Max;

				throw std::runtime_error("Bereichsueberschreitung in Histo1D!!!"+ temp.str());
			}

			if ((x1 < Min))
			{
				throw std::runtime_error("Bereichsunterschreitung in Histo1D!!!");
			}

			Werte[(int) floor((wert-Min)/dI)]++;

			if(MaximumValue < Werte[(int) floor((wert-Min)/dI)])
				MaximumValue = Werte[(int) floor((wert-Min)/dI)];
		}

		double GetNrInBin(int bin)
		{
			return Werte[bin];
		}

		double GetNrInBinNormiert(int bin)
		{
			return Werte[bin]/N;
		}

		double GetNrInBinNormiertMaximum(int bin)
		{
			return Werte[bin]/MaximumValue;
		}

		double GetRangeInBin(int bin)
		{
			return (Min+(bin)*dI + dI/2.0);
		}

		int GetNrBins()
		{
			return NrBins;
		}

		double GetCumulativeNrInBin(int bin)
			{
				double cumValue=0.0;

				for(int i=0; i<=bin;i++)
					cumValue+=Werte[i];

				return cumValue;
			}

		double GetCumulativeNrInBinNormiert(int bin)
			{
				double cumValue=0.0;

				for(int i=0; i<=bin;i++)
					cumValue+=Werte[i];

				return cumValue/N;
			}


		double GetIntervallThickness()
		{
			return dI;
		}

		int GetNrCounts()
		{
			return N;
		}

		void ResetHistogram1D(double min, double max, int nrbins, double initValue)
		{
			delete [] Werte;

			NrBins = nrbins;


					N= 0;
					Min = min;


					Max = max;


					wert = 0;


					dI = (max-min)/(1.0*nrbins);



					MaximumValue = 0;
					//Werte = new double[NrBins1+1][NrBins2+1];

					Werte = new double[NrBins+1];
					for(int i = 0; i < (NrBins+1); ++i) {
						Werte[i] = initValue;
					}
		}

		void operator += (Histogram1D& rhsHG)
			{
				for(int bin = 0; bin < rhsHG.GetNrBins() ;bin++)
				  	{
						Werte[bin] += rhsHG.GetNrInBin(bin);

				  	}

				N += rhsHG.GetNrCounts();
			}

		const double GetNrAtValue(double _value) const
				{
					return Werte[(int) floor((wert-Min)/dI)];
				}

private:


	int NrBins; //Anzahl der Bins
	int N; //Anzahl aller Eintraege
	double Min; // untere Grenze
	double Max; // obere Grenze
	double wert; //Wert weiterreichen
	double dI; //Intervalleinteilung

	double* Werte; //Speicherung im Intervall

	double MaximumValue; //Maximalwert der Stichprobe

};

#endif //Histogram1D.h
