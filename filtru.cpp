#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>

using namespace std;

int main(int argc, char * argv[]){

	int rank;
	int n_proc;
	string line;
	int i, j;
	int lines_proc = 0;


	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);

	// Citim argumentele si le separam  

	char argv1[20], argv2[20],argv3[20];
	char *w = strtok (argv[6]," ");

	if(w != NULL){
		strcpy (argv1, w + 10);
		w = strtok(NULL, " ");
		if (w != NULL){
			strcpy (argv2, w);
			w = strtok(NULL, " ");
			if (w != NULL){
				strcpy (argv3, w);
			}
		}
	}

	// Citim topologia grafului
	ifstream file(argv1);

	// ne intereseaza doar vecinii nodului rank
	for (i = 0; i <= rank; i++){
		getline (file, line);
	}
	
	// Extragem vecinii din linia citita
	char *words = new char [line.length()+1];
	strcpy (words, line.c_str());
	
	char *p = strtok (words," ");
	p = strtok(NULL, " ");

	vector<int> neighbors;
	char x[10];
	while (p != NULL){
		strcpy(x, p);
		neighbors.push_back(atoi(x));
		p = strtok(NULL, " ");
	}
	delete[] words;

	
	// Alg heartbeat

	vector<int> children;
	int parent = -1;
	int element;


	// Plecam din radacina si trimitem vecinilor rank-ul si asteptam 
	// confirmarea ca rank-ul curent este parintele lor
	if (rank == 0){
		for (i = 0; i < neighbors.size(); i++){

			MPI_Send( &rank, 1, MPI_INT, 
					neighbors[i], 1, MPI_COMM_WORLD);

			MPI_Recv(&element, 1, MPI_INT,
					MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);

			if (element == rank){
				children.push_back(status.MPI_SOURCE);
			}
		}
	}

	// Trimit rank-ul si astept o confirmare. Daca are alt parinte, atunci
	// il sterg din lista de copii
	else{
		// Astept ca un proces sa-mi spuna ca este parintele meu si trimit
		// inapoi ca raspuns parintele meu
		MPI_Recv(&element, 1, MPI_INT,MPI_ANY_SOURCE,
				1, MPI_COMM_WORLD, &status);

		if (parent == -1){
			parent = element;
		}
		
		MPI_Send(&parent, 1, MPI_INT, status.MPI_SOURCE,
				 1, MPI_COMM_WORLD);

		// Trimit rank-ul si astept o confirmare. Daca copilul raspunde cu
		// rank-ul meu, atunci il adaug in lista de copii
		for (i = 0; i < neighbors.size(); i++){

			if (neighbors[i] != parent){
				MPI_Send(&rank, 1, MPI_INT, neighbors[i],
						1, MPI_COMM_WORLD);

				MPI_Recv(&element, 1, MPI_INT, MPI_ANY_SOURCE, 1,
					MPI_COMM_WORLD, &status);

				if (element == rank){
					children.push_back(status.MPI_SOURCE);
				}
			}
		}
	}

	// stergem vectorul de vecini
	neighbors.clear();

	MPI_Barrier(MPI_COMM_WORLD);

	// Deschidem fisierul imagini.in
	ifstream images(argv2);
	int num_img, k;
	int filtru;
	string str, img_in, img_out;
	int n, m;
	int** matrix;
	int div;
	int rest;

	images >> num_img;

	// Citim lista imaginilor de procesat

	for (k = 0; k < num_img; k++){

		images >> str >> img_in >> img_out;

		// Setam filtrul 
		if (str.compare("smooth") == 0){
			filtru = 1;
		}
		else if (str.compare("blur") == 0){
			filtru = 2;
		}
		else if (str.compare("sharpen") == 0){
			filtru = 3;
		}
		else if (str.compare("mean_removal") == 0){
			filtru = 4;
		}

		ifstream in(img_in.c_str());
		ofstream out(img_out.c_str());

		// Radacina va citi matricea imaginii si o va trimite mai departe
		if (rank == 0){

			getline (in, line); out << line << '\n';
			getline (in, line); out << line << '\n';

			in >> m >> n;
			out << m << ' ' << n << '\n';

			in >> str; out << str << '\n';

			// adaugam bordarea
			n += 2;
			m += 2;

			// matricea imaginii
			matrix = new int* [n];
			for(i = 0; i < n; i++) { 
			   matrix[i] = new int[m];
			}

			// bordam cu 0 
			for(i = 0; i < n; i++) { 
			   matrix[i][0] = 0;
			   matrix[i][m-1] = 0;
			}
			for(i = 0; i < m; i++) { 
			   matrix[0][i] = 0;
			   matrix[n-1][i] = 0;
			}

			// citim matricea imaginii
			for (i = 1; i < n-1; i++){
				for (j = 1; j < m-1; j++){
					in >> matrix[i][j];
				}
			}


			// calculam numarul de fasii ce vor fi trimise
			if (children.size() > 0){
				// Avem n-2 linii ce trebuiesc procesate 
				div = (n-2) / children.size();
				div += 2;
				rest = (n-2) % children.size();
			}

			int index = 0;

			// nr de linii este mai mic decat cel al copiilor
			// trimit cate o linie copiilor pana le epuizez, iar restul
			// copiilor le trimit -1 pentru a-i anunta ca nu primesc nimic
			if (n-2 < children.size()  && children.size() > 0){
				div = 3;

				for (i = 0; i < n-2; i++){

					MPI_Send(&div, 1, MPI_INT, children[i], 
							1, MPI_COMM_WORLD);

					MPI_Send(&m, 1, MPI_INT, children[i], 
							1, MPI_COMM_WORLD);

					for (j = 0; j < div; j++){

						MPI_Send(matrix[index + j], m, MPI_INT, 
							children[i], 1, MPI_COMM_WORLD);

					}

					index++;
				}

				for (i = n-2; i < children.size(); i++){
					div = -1;
					MPI_Send(&div, 1, MPI_INT, children[i],
							 1, MPI_COMM_WORLD);
				}			
			}

			// nr de linii este mai mare decat cel al copiilor
			else
			for (i = 0; i < children.size(); i++){
				
				// ultimul copil primeste si restul de linii
				// ramase in urma divizarii
				if (i == children.size() - 1){
					div += rest;

				}

				// trimit copiilor nr de linii ale fasiei, m si fasia 
				MPI_Send(&div, 1, MPI_INT, children[i], 
						children[i], MPI_COMM_WORLD);
				MPI_Send(&m, 1, MPI_INT, children[i],
						children[i], MPI_COMM_WORLD);

				for (j = 0; j < div; j++){

					MPI_Send(matrix[index + j], m, MPI_INT, 
						     children[i], children[i], MPI_COMM_WORLD);

				}

				index = index + div - 2;
			
			}			
		}
		//celelalte noduri
		else {
			// primesc numarul de linii
			MPI_Recv(&n, 1, MPI_INT, parent, rank,
					 MPI_COMM_WORLD, &status);

			// Daca am primit -1, inseamna ca nu voi mai primi nicio fasie,
			// asa ca anunt copiii ca nu primesc nici ei
			if(n == -1){
				for (i = 0; i < children.size(); i++){
					MPI_Send(&n, 1, MPI_INT, children[i],
							 1, MPI_COMM_WORLD);
				}
			}
			// Primesc fasia
			else{

			MPI_Recv(&m, 1, MPI_INT, parent, 
			         rank, MPI_COMM_WORLD, &status);
			
			// alocam matricea pentru noua fasie
			matrix = new int* [n];
			for(i = 0; i < n; i++) { 
			   matrix[i] = new int [m];
			}

			// primesc matricea
			for (i = 0; i < n; i++){
				MPI_Recv(matrix[i], m, MPI_INT, parent, rank,
						 MPI_COMM_WORLD, &status);
			}

			// impartim matricea la numarul de copii
			if (children.size() > 0){
				div = (n-2) / children.size();
				div += 2;
				rest = (n-2) % children.size();
			}

			int index = 0;

			// nr de linii este mai mic decat cel al copiilor.
			// trimit cate o linie copiilor pana le epuizez, iar restul
			// copiilor le trimit -1 pentru a-i anunta ca nu primesc nimic
			if (n-2 < children.size()  && children.size() > 0){
				div = 3;

				for (i = 0; i < n-2; i++){

					MPI_Send(&div, 1, MPI_INT, children[i],
							1, MPI_COMM_WORLD);

					MPI_Send(&m, 1, MPI_INT, children[i],
							 1, MPI_COMM_WORLD);

					for (j = 0; j < div; j++){
						MPI_Send(matrix[index + j], m, MPI_INT, 
							children[i], 1, MPI_COMM_WORLD);

					}

					index++;
				}

				for (i = n-2; i < children.size(); i++){
					div = -1;
					MPI_Send(&div, 1, MPI_INT, children[i], 
							1, MPI_COMM_WORLD);
				}
			}

			// nr de linii este mai mare decat cel al copiilor
			else
			for (i = 0; i < children.size(); i++){

				// ultimul copil primeste si restul de linii
				// ramase in urma divizarii
				if (i == children.size() - 1){
					div += rest;
				}

				// trimit copiilor fasii egale din matricea primita
				MPI_Send(&div, 1, MPI_INT, children[i],
				 		children[i], MPI_COMM_WORLD);
				MPI_Send(&m, 1, MPI_INT, children[i], 
						children[i], MPI_COMM_WORLD);

				for (j = 0; j < div; j++){
					MPI_Send(matrix[index + j], m, MPI_INT,
						children[i], children[i], MPI_COMM_WORLD);
				}

				index = index + div - 2;
			}
		}

		} // end else (alt rank)*/



		// ******** APLICAM FILTRUL  **********


		// verificam daca suntem pe o frunza
		if (children.size() == 0 && n > 0){

			int fa,fb,fc,fd,fe,ff,fg,fh,fi;
			// coeficientii filtrelor

			if(filtru == 1){
				fa = 1;
				fb = 1;
				fc = 1;
				fd = 1;
				fe = 1;
				ff = 1;
				fg = 1;
				fh = 1;
				fi = 1;
			}
			else if (filtru == 2){
				fa = 1;
				fb = 2;
				fc = 1;
				fd = 2;
				fe = 4;
				ff = 2;
				fg = 1;
				fh = 2;
				fi = 1;
			}
			else if(filtru == 3){
				fa = 0;
				fb = -2;
				fc = 0;
				fd = -2;
				fe = 11;
				ff = -2;
				fg = 0;
				fh = -2;
				fi = 0;
			}
			else if(filtru == 4){
				fa = -1;
				fb = -1;
				fc = -1;
				fd = -1;
				fe = 9;
				ff = -1;
				fg = -1;
				fh = -1;
				fi = -1;
			}

			// Noua matrice a imaginii
			int ** result;

			result = new int* [n];
			for(i = 0; i < n; i++) { 
			   result[i] = new int [m];
			}

			// Aplicam filtrul pe fiecare element din matricea bordata
			for (i = 1; i < n-1; i++){
				for (j = 1; j < m-1; j++){
					result[i][j] = fa * matrix[i-1][j-1];
					result[i][j] += fb * matrix[i-1][j];
					result[i][j] += fc * matrix[i-1][j+1];
					result[i][j] += fd * matrix[i][j-1];
					result[i][j] += fe * matrix[i][j];
					result[i][j] += ff * matrix[i][j+1];
					result[i][j] += fg * matrix[i+1][j-1];
					result[i][j] += fh * matrix[i+1][j];
					result[i][j] += fi * matrix[i+1][j+1];

					if (filtru == 1){
						result[i][j] /= 9;
					}
					else if (filtru == 2){
						result[i][j] /= 16;
					}
					else if (filtru == 3){
						result[i][j] /= 3;
					}
					else if (filtru == 4){
						result[i][j] /= 1;
					}

					if (result[i][j] < 0){
						result[i][j] = 0;
					}
					else if (result[i][j] > 255){
						result[i][j] = 255;
					}
				}
			}

			// numarul de linii procesate
			lines_proc += (n-2);		

			// trimit parintelui numarul de linii
			MPI_Send(&n, 1, MPI_INT, parent, rank, MPI_COMM_WORLD);

			// trimit noua matrice incepand cu linia 1
			for (i = 1; i < n; i++){
				MPI_Send(result[i], m, MPI_INT, 
						parent, rank, MPI_COMM_WORLD);
			}

			// DEALOCARE matrice de rezultate
			for (i = 0; i < n; i++){
				delete[] result[i];
			}
			delete[] result;
			
		}

		// procesul este nod intermediar => construieste noua matrice
		else{
			j = 1;

			for (i = 0; i < children.size(); i++){
				int num;
				// primesc nr de linii si stim ca vor fi trimise n-1 linii
				// incepand cu a doua
				
				MPI_Recv(&num, 1, MPI_INT, children[i],
						children[i], MPI_COMM_WORLD, &status);
				
				num --;

				// Actualizam matricea mare
				for (int ii = 0; ii < num; ii++){
					MPI_Recv(matrix[j], m, MPI_INT, children[i],
							children[i], MPI_COMM_WORLD, &status);
					j++;
				}

				// Scadem j deoarece in matrice ultima linie
				// adaugata face parte din bordarea fasiei
				j--;
			}

			// Daca nu s-a ajuns in nodul radacina, atunci trimitem
			// matricea(fasia) mai departe la parinte
			if (parent != -1){

				MPI_Send(&n, 1, MPI_INT, parent, rank, MPI_COMM_WORLD);

				for (i = 1; i < n; i++){
					MPI_Send(matrix[i], m, MPI_INT, parent,
							rank, MPI_COMM_WORLD);
				}
			}
		}

		// Asteptam ca toate procesele sa termine lucrul cu matricea imaginii
		MPI_Barrier(MPI_COMM_WORLD);

		// scriem in fisier noua imagine rezultata in urma filtrului
		if (rank == 0){
			for (i = 1; i < n-1; i++){
				for (j = 1; j < m-1; j++){
					out << matrix[i][j] << '\n';
				}
			}
		}
		out.close();
		in.close();

	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Vectori pentru statistica: primul este acumulator, iar in al doilea
	// primesc statistici de la vecini 
	int* statistics = new int [n_proc];
	int* Rstatistics = new int [n_proc];
	
	// Initilizam vectorii
	for (i = 0 ; i < n_proc; i++){
		statistics[i] = 0;
		Rstatistics[i] = 0;
	}

	// Procesul radacina va trimite mesaj catre copii, si va
	// astepta vectorii de statistici ai acestora
	if (rank == 0){
		for (i = 0; i < children.size(); i++){
			MPI_Send(statistics, n_proc, MPI_INT, children[i],
					 rank, MPI_COMM_WORLD);
		}

		// primim statistici de la vecini si actualizam
		// vectorul de statistici al procesului curent
		for (i = 0; i < children.size(); i++){
			MPI_Recv(Rstatistics, n_proc, MPI_INT,
					MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);

			for (j = 0; j < n_proc; j++){
				statistics[j] += Rstatistics[j];
			}
		}
	}

	// Un alt proces
	else{

		// Proces frunza
		if (children.size() == 0){
			// primesc mesaj de la parinte
			MPI_Recv(Rstatistics, n_proc, MPI_INT, parent,
					 parent, MPI_COMM_WORLD, &status);

			// trimit inapoi statistica cu nr de linii procesate
			statistics[rank] = lines_proc;
			MPI_Send(statistics, n_proc, MPI_INT, parent, 
					parent, MPI_COMM_WORLD);
		}

		// nod intermediar
		else{
			// primesc mesaj de la parinte
			MPI_Recv(Rstatistics, n_proc, MPI_INT, parent, parent, 
					MPI_COMM_WORLD, &status);

			// il trimit la toti copiii
			for (i = 0; i < children.size(); i++){
				MPI_Send(Rstatistics, n_proc, MPI_INT, children[i],
						rank, MPI_COMM_WORLD);
			}

			// primesc de la copii statisticile si le cumulez
			for (i = 0; i < children.size(); i++){
				MPI_Recv(Rstatistics, n_proc, MPI_INT,
					MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);

				for (j = 0; j < n_proc; j++){
					statistics[j] += Rstatistics[j];
				}
			}

			// trimit statisticile cumulate parintelui
			MPI_Send(statistics, n_proc, MPI_INT, parent,
					 parent, MPI_COMM_WORLD);
		}		
	}

	// Radacina va scrie statisticile primite
	if (rank == 0){
		ofstream stat_file(argv3);
		for (i = 0 ; i < n_proc; i++){
			stat_file << i << ": " << statistics[i] <<'\n';
		}
		stat_file.close();
	}

	// DEALOCARE matrice imagine
	for (i = 0; i < n; i++){
		delete[] matrix[i];
	}
	delete[] matrix;


	
	images.close();
	file.close();

	MPI_Finalize();
	
	return 0;
}
