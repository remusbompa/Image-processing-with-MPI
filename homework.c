#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

typedef struct {
	int type;
	int width;
	int height;
	unsigned char maxval;
	void *pixels;
}image;

typedef struct{
	unsigned char r,g,b;
}__attribute__((packed)) ColorP;

typedef unsigned char GrayscaleP;

typedef struct{
	int start;
	int end;
}GroupP;

image in, out;
int nrFilters;
char** nameFilters;
int rank;
int nProcesses;

void readInput(const char * fileName, image *img) {
	FILE* fin = fopen( fileName, "r");
	fscanf(fin, "P%d\n", &img -> type);
	fscanf(fin, "%d %d\n", &img -> width, &img -> height);
	fscanf(fin, "%hhu\n", &img -> maxval);
	if(img -> type == 6)
	{
		img -> pixels = malloc( sizeof(ColorP*) * img -> height);
		int i;
		for( i=0; i< img -> height; i++ ){
			((ColorP**)img -> pixels)[i] = (ColorP*)malloc( sizeof(ColorP) * img -> width);
			int j;
			for( j=0; j< img -> width; j++){
				ColorP* pixel= &((ColorP**)img -> pixels)[i][j];
				unsigned char culori[3];
				fread( culori, 1, 3, fin);
				pixel -> r = culori[0]; 
				pixel -> g = culori[1];
				pixel -> b = culori[2];
			}
		}
	}
	else if(img -> type == 5)
	{
		img -> pixels = malloc( sizeof(GrayscaleP*) * img -> height);
		int i;
		for( i=0; i< img -> height; i++ ){
			((GrayscaleP**)img -> pixels)[i] = (GrayscaleP*)malloc( sizeof(GrayscaleP) * img -> width);
			int j;
			for( j=0; j< img -> width; j++){
				GrayscaleP* pixel= &((GrayscaleP**)img -> pixels)[i][j];
				fread( pixel, 1, 1, fin);
			}
		}
	}
}

void writeData(const char * fileName, image *img) {

	FILE* fout = fopen( fileName, "wb");
	
	fprintf(fout, "P%d\n", img -> type);
	fprintf(fout, "%d %d\n", img -> width, img -> height);
	fprintf(fout, "%d\n", img -> maxval);
	int i,j;

	if( img -> type == 6)
	{
		for(i= 0; i< img -> height; i++){
			for( j=0; j< img -> width; j++){
				ColorP* pixel= &((ColorP**)img -> pixels)[i][j];
				unsigned char culori[3];
				culori[0] = pixel -> r; 
				culori[1] = pixel -> g; 
				culori[2] = pixel -> b; 
				fwrite( culori, 1, 3, fout);
			}
		}
	}
	else if( img -> type == 5){
		for(i= 0; i< img -> height; i++){
			for( j=0; j< img -> width; j++){
				GrayscaleP* pixel= &((GrayscaleP**)img -> pixels)[i][j];
				fwrite( pixel, 1, 1, fout);
			}
		}
	}
}

void aplicaFiltruColor(ColorP* newPixel, float k[3][3], ColorP vecini[3][3]){
	float aux_r = 0, aux_g = 0, aux_b = 0;
	int i,j;
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			float vecin_r = vecini[i][j].r;
			vecin_r = vecin_r * k[i][j];
			float vecin_g = vecini[i][j].g;
			vecin_g = vecin_g * k[i][j];
			float vecin_b = vecini[i][j].b;
			vecin_b = vecin_b * k[i][j];

			aux_r += vecin_r;
			aux_g += vecin_g;
			aux_b += vecin_b;

		}
	}
	
	newPixel->r = (unsigned char)aux_r;
	newPixel->g = (unsigned char)aux_g;
	newPixel->b = (unsigned char)aux_b;
}

void aplicaFiltruGrayscale(GrayscaleP* newPixel, float k[3][3], GrayscaleP vecini[3][3]){
	float aux = 0;
	int i,j;
	for(i = 0; i < 3; i++){
		for(j = 0; j < 3; j++){
			float vecin = vecini[i][j];
			vecin = vecin * k[i][j];
			aux += vecin;
		}
	}
	*newPixel = (unsigned char)aux;
}

void copyMatrix(int tip, void* pixels, void* buf, int start, int end, int width){
	int istart = start / width, jstart = start % width;
	int iend = end / width, jend = end % width;
	if(tip == 6){
		ColorP** colorPixels = (ColorP**)pixels;
		ColorP* buffer = (ColorP*)buf;
		int i=istart, j=jstart;
		while(1)
		{
			int ind = i * width + j - start;
			buffer[ind].r = colorPixels[i][j].r;
			buffer[ind].g = colorPixels[i][j].g;
			buffer[ind].b = colorPixels[i][j].b;
			if( i == iend && j == jend) break;
			j++;
			if(j == width)
			{
				j=0;
				i++;
			}			
		}
	}else{
		GrayscaleP** grayscalePixels = (GrayscaleP**)pixels;
		GrayscaleP* buffer = (GrayscaleP*)buf;
		int i=istart, j=jstart;
		while(1)
		{
			int ind = i * width + j - start;
			buffer[ind] = grayscalePixels[i][j];
			if( i == iend && j == jend) break;
			j++;
			if(j == width)
			{
				j=0;
				i++;
			}			
		}
	}
}

void copyBuffer(int tip, void* pixels, void* buf, int start, int end, int width){
	int istart = start / width, jstart = start % width;
	int iend = end / width, jend = end % width;
	if(tip == 6){
		ColorP** colorPixels = (ColorP**)pixels;
		ColorP* buffer = (ColorP*)buf;
		int i=istart, j=jstart;
		while(1)
		{
			int ind = i * width + j - start;
			colorPixels[i][j].r = buffer[ind].r;
			colorPixels[i][j].g = buffer[ind].g;
			colorPixels[i][j].b = buffer[ind].b;
			if( i == iend && j == jend) break;
			j++;
			if(j == width)
			{
				j=0;
				i++;
			}			
		}
	}else{
		GrayscaleP** grayscalePixels = (GrayscaleP**)pixels;
		GrayscaleP* buffer = (GrayscaleP*)buf;
		int i=istart, j=jstart;
		while(1)
		{
			int ind = i * width + j - start;
			grayscalePixels[i][j] = buffer[ind];
			if( i == iend && j == jend) break;
			j++;
			if(j == width)
			{
				j=0;
				i++;
			}			
		}
	}
}

void computePixels(int* param, int rank){
	MPI_Status status;

	int start = param[0];
	int end = param[1];
	int tip = param[2];
	int width = param[3];
	int height = param[4];

	float k[5][3][3] ={
				{{1.0f/9.0f,1.0f/9.0f,1.0f/9.0f},{1.0f/9.0f,1.0f/9.0f,1.0f/9.0f},{1.0f/9.0f,1.0f/9.0f,1.0f/9.0f}}, 
				{{1.0f/16.0f,2.0f/16.0f,1.0f/16.0f},{2.0f/16.0f,4.0f/16.0f,2.0f/16.0f},{1.0f/16.0f,2.0f/16.0f,1.0f/16.0f}},
				{{0,-2.0f/3.0f,0},{-2.0f/3.0f,11.0f/3.0f,-2.0f/3.0f},{0,-2.0f/3.0f,0}},
				{{-1,-1,-1},{-1,9,-1},{-1,-1,-1}},
				{{0,1,0},{0,0,0},{0,-1,0}} 
				};


	int istart = start / width, jstart = start % width;
	int iend = end / width, jend = end % width;

	int sentStart = MAX(start - width - 1, 0);
	int sentEnd = MIN(end + width + 1, width * height - 1);

	int nr = sentEnd - sentStart + 1;
	if(tip == 6){
		ColorP *oldPixels = (ColorP*)malloc(nr * sizeof(ColorP));


		//primeste pixelii de la sentStart la sentEnd in oldPixels
		if(rank == 0) copyMatrix(in.type, in.pixels, oldPixels, sentStart, sentEnd, width); 
		else MPI_Recv(oldPixels, 3 * nr, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);


		//pentru fiecare filtru
		int x,f;
		for(x = 0; x < nrFilters; x++){
			if(!strcmp(nameFilters[x], "smooth")) f =0;
			else if(!strcmp(nameFilters[x], "blur")) f =1;
			else if(!strcmp(nameFilters[x], "sharpen")) f =2;
			else if(!strcmp(nameFilters[x], "mean")) f =3;
			else if(!strcmp(nameFilters[x], "emboss")) f =4;


			//construire noi pixeli pentru filtrul dat
			ColorP* newPixels = (ColorP*)malloc(nr * sizeof(ColorP));
			//parcurgere pixeli de la start la end (inclusiv)
			int i=istart, j=jstart;
			while(1)
			{
				//indicele celulei mele in newPixels si oldPixels
				int ind = i * width + j - sentStart;
				//daca pixelul e pe margine ramane neschimbat
				if(i == 0 || j == 0 || i == height - 1 || j == width - 1){
					newPixels[ind].r = oldPixels[ind].r;
					newPixels[ind].g = oldPixels[ind].g;
					newPixels[ind].b = oldPixels[ind].b;
				}
				//altfel se aplica filtrul folosind matricea vecinilor
				else{
					ColorP vecini [3][3] = { {oldPixels[ind - width - 1], oldPixels[ind - width], oldPixels[ind - width + 1]} ,
												    {oldPixels[ind - 1], oldPixels[ind], oldPixels[ind + 1]},
											 		{oldPixels[ind + width - 1], oldPixels[ind + width], oldPixels[ind + width + 1]}, 
											};
					aplicaFiltruColor(&newPixels[ind], k[f], vecini);
				}


				if( i == iend && j == jend) break;
				j++;
				if(j == width)
				{
					j=0;
					i++;
				}
			}

			//primesc mesaj din stanga
			if(rank >= 1){
				MPI_Recv(&newPixels[0], 3 * (width + 1), MPI_CHAR, rank - 1, x+1, MPI_COMM_WORLD, &status);	
			}
			//trimit mesaj in dreapta
			if(rank < nProcesses - 1){
				MPI_Send(&newPixels[end - sentStart - width ], 3 * (width +1), MPI_CHAR, rank + 1, x+1, MPI_COMM_WORLD);
			}
			//primesc mesaj din dreapta
			if(rank < nProcesses - 1){
				MPI_Recv(&newPixels[end -sentStart +1 ], 3 * (width + 1), MPI_CHAR, rank + 1, x+1, MPI_COMM_WORLD, &status);	
			}
			//trimit mesaj in stanga	
			if(rank >= 1){
				MPI_Send(&newPixels[start - sentStart], 3 * (width + 1), MPI_CHAR, rank - 1, x+1, MPI_COMM_WORLD);
			}

			free(oldPixels);
			oldPixels = newPixels;	

		}

	//dupa ultimul filtru, procesele trimit pixelii de la start la end procesati spre procesul 0
		if(rank == 0){
			copyBuffer(in.type, out.pixels, &oldPixels[start - sentStart], start, end, width); 
		}else{
			MPI_Send(&oldPixels[start - sentStart], 3 * (end - start + 1), MPI_CHAR, 0, nrFilters + 1, MPI_COMM_WORLD);
		}


	}else{
	
		GrayscaleP* oldPixels = (GrayscaleP*)malloc(nr * sizeof(GrayscaleP));
		if(rank == 0) copyMatrix(in.type, in.pixels, oldPixels, sentStart, sentEnd, width); 
		else MPI_Recv(oldPixels, nr, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
		//pentru fiecare filtru
		int x,f;
		for(x = 0; x < nrFilters; x++){
			if(!strcmp(nameFilters[x], "smooth")) f =0;
			else if(!strcmp(nameFilters[x], "blur")) f =1;
			else if(!strcmp(nameFilters[x], "sharpen")) f =2;
			else if(!strcmp(nameFilters[x], "mean")) f =3;
			else if(!strcmp(nameFilters[x], "emboss")) f =4;


			//construire noi pixli pentru filtrul dat
			GrayscaleP* newPixels = (GrayscaleP*)malloc(nr * sizeof(GrayscaleP));
			//parcurgere pixeli de la start la end (inclusiv)
			int i=istart, j=jstart;
			while(1)
			{
				int ind = i * width + j - sentStart;
				//daca pixelul e pe margine ramane neschimbat
				if(i == 0 || j == 0 || i == height - 1 || j == width - 1){
					newPixels[ind] = oldPixels[ind];
				}
				//altfel se aplica filtrul folosind matricea vecinilor
				else{
					GrayscaleP vecini [3][3]= {{oldPixels[ind - width - 1], oldPixels[ind - width],oldPixels[ind - width + 1]},
														{oldPixels[ind - 1], oldPixels[ind], oldPixels[ind + 1]},
													{oldPixels[ind + width - 1], oldPixels[ind + width], oldPixels[ind + width + 1]}
													};
					aplicaFiltruGrayscale(&newPixels[ind], k[f], vecini);
				}


				if( i == iend && j == jend) break;
				j++;
				if(j == width)
				{
					j=0;
					i++;
				}
			}
			
			//primesc mesaj din stanga
			if(rank >= 1){
				MPI_Recv(&newPixels[0], width + 1, MPI_CHAR, rank - 1, x+1, MPI_COMM_WORLD, &status);	
			}
			//trimit mesaj in dreapta
			if(rank < nProcesses - 1){
				MPI_Send(&newPixels[end - sentStart - width ], width +1, MPI_CHAR, rank + 1, x+1, MPI_COMM_WORLD);
			}
			//primesc mesaj din dreapta
			if(rank < nProcesses - 1){
				MPI_Recv(&newPixels[end -sentStart +1 ], width + 1, MPI_CHAR, rank + 1, x+1, MPI_COMM_WORLD, &status);	
			}
			//trimit mesaj in stanga	
			if(rank >= 1){
				MPI_Send(&newPixels[start - sentStart], width + 1, MPI_CHAR, rank - 1, x+1, MPI_COMM_WORLD);
			}
			free(oldPixels);
			oldPixels = newPixels;		
		}

	//dupa ultimul filtru, procesele trimit pixelii de la start la end procesati spre procesul 0
		if(rank == 0){
			copyBuffer(in.type, out.pixels, &oldPixels[start- sentStart], start, end, width);
		}else{
			MPI_Send(&oldPixels[start - sentStart], (end - start + 1), MPI_CHAR, 0, nrFilters + 1, MPI_COMM_WORLD);
		}

	}
}

int main(int argc, char * argv[]) {

	MPI_Init(&argc, &argv);

	nrFilters = argc - 3;
	nameFilters = &argv[3];

	if(!strcmp(nameFilters[0], "bssembssem")){
		nrFilters = 10;
		nameFilters = malloc(nrFilters * sizeof(char*));
		int i;
		for(i = 0; i < nrFilters; i++){
			nameFilters[i] = malloc(20 * sizeof(char));
		}
		strcpy(nameFilters[0], "blur");
		strcpy(nameFilters[1], "smooth");
		strcpy(nameFilters[2], "sharpen");
		strcpy(nameFilters[3], "emboss");
		strcpy(nameFilters[4], "mean");
		strcpy(nameFilters[5], "blur");
		strcpy(nameFilters[6], "smooth");
		strcpy(nameFilters[7], "sharpen");
		strcpy(nameFilters[8], "emboss");
		strcpy(nameFilters[9], "mean");
	}

	MPI_Status status;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

	if(rank == 0){
		readInput(argv[1], &in);

		int height = in.height;
 		int width = in.width;
		int nr_grupuriP = height * width;

		//impartire pixeli din matrice in P procese
		GroupP group_proc [nProcesses];
		int chunks = nr_grupuriP / nProcesses;
		int reminder = nr_grupuriP % nProcesses;

		int p = 0;
		int i=0;
		while( reminder > 0 )
		{
			GroupP* agrup = &group_proc[i];
			agrup -> start = p;
			agrup -> end = p + chunks - 1 + 1;
			reminder--;
			p = agrup -> end + 1;
			i++;

		}
		while( i < nProcesses ){
			GroupP* agrup = &group_proc[i];
			agrup -> start = p;
			agrup -> end = p + chunks - 1;
			p = agrup -> end + 1;
			i++;
		}

		for(i =1; i < nProcesses; i++){
			GroupP* agrup = &group_proc[i];
			int start = agrup->start, end = agrup->end;
			int tip = in.type;
			int param [5] = {start, end, tip, width, height};

			MPI_Send(param, 5, MPI_INT, i, 1, MPI_COMM_WORLD);


			int sentEnd = MIN( (end + width + 1), height * width - 1);
			int sentStart = MAX( (start - width - 1), 0);
			int nr = sentEnd - sentStart + 1;

			

			if(tip == 6){
				ColorP* buffer = (ColorP*)malloc(nr * sizeof(ColorP));
				copyMatrix(in.type, in.pixels, buffer, sentStart, sentEnd, width);
				MPI_Send(buffer, 3 * nr, MPI_CHAR, i, 0, MPI_COMM_WORLD);
				free(buffer);
			}
			else{
				GrayscaleP* buffer = (GrayscaleP*)malloc(nr * sizeof(GrayscaleP));
				copyMatrix(in.type, in.pixels, buffer, sentStart, sentEnd, width);
				MPI_Send(buffer, nr, MPI_CHAR, i, 0, MPI_COMM_WORLD);
				free(buffer);
			}



		}

		//procesul 0 calculeaza primul grup de pixeli
		GroupP* agrup = &group_proc[0];
		int start = agrup->start, end = agrup->end;
		int tip = in.type;
		int param [5] = {start, end, tip, width, height};
		
		//construire imagine finala
		out.height = height;
		out.width = width;
		out.type = in.type;
		out.maxval = in.maxval;
		if( in.type == 6){
			out.pixels = malloc( height * sizeof(ColorP*));
			ColorP** pixels = ((ColorP**)out.pixels);
			for(  i = 0; i < height; i++){
				pixels[i] = (ColorP*) malloc( width * sizeof(ColorP));
			}
		}
		else if( in.type == 5){
			out.pixels = malloc( height * sizeof(GrayscaleP*));
			GrayscaleP** pixels = ((GrayscaleP**)out.pixels);
			for(  i = 0; i < height; i++){
				pixels[i] = (GrayscaleP*) malloc( width * sizeof(GrayscaleP));
			}
		}

		computePixels(param, rank);

		//primire pixeli procesati de celelalte procese

		for(p = 1; p < nProcesses; p++){
			GroupP* agrup = &group_proc[p];
			int start = agrup->start, end = agrup->end;

			if(in.type == 6){
				ColorP* buffer = (ColorP*)malloc((end - start + 1) * sizeof(ColorP));
				MPI_Recv(buffer, 3 * (end - start + 1), MPI_CHAR, p, nrFilters + 1, MPI_COMM_WORLD, &status);
				copyBuffer(out.type, out.pixels, buffer, start, end, width);
				free(buffer);
			}else{
				GrayscaleP* buffer = (GrayscaleP*)malloc((end - start + 1) * sizeof(GrayscaleP));
				MPI_Recv(buffer, (end - start + 1), MPI_CHAR, p, nrFilters + 1, MPI_COMM_WORLD, &status);
				copyBuffer(out.type, out.pixels, buffer, start, end, width);
				free(buffer);
			}
		}
		//scriere fisier de iesire
		writeData(argv[2], &out);
	}else{
		int param[5];
		MPI_Recv(param, 5, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		computePixels(param, rank);
	}

	MPI_Finalize();
	return 0;
}
