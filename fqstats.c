// gcc fqstats.c -o fqstats -lm

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double LOOKUPTABLE[100];


typedef struct 
{
	char header[1024];
	char dna[2097152];
	char third[1024];
	char qua[2097152];
}record;


double recordStats(record *rec)
{
	int pos=0;
	double sum_prob=0;
	int q;
	int i=0;
	int gc=0;

	double meanqscore=0;

	while (q=rec->qua[pos++])
	{
		if (q=='\n')
			break;
		q-=33;
		sum_prob+=LOOKUPTABLE[q];
	}

	pos--;

	for (i=0;rec->dna[i]!='\n';i++)
		switch(rec->dna[i])
		{
			case 'c':
			case 'C':
			case 'g':
			case 'G':
				gc++;
				break;
		}


	if (!sum_prob) meanqscore=0;
	else{

		double fmean = sum_prob/(pos);
		meanqscore = -10 * log10(fmean);
	}

	for (i=0;rec->header[i]!='\0';i++)
	{
		if (rec->header[i]==' ')
		{
			rec->header[i]='\0';
			break;
		}
	}

	printf("%s\t%d\t%f\t%f\n",(rec->header)+1,pos,meanqscore,((double)gc)/pos);
}

int main (int argc, char **argv)
{

	unsigned long int pos=0;

	FILE * fp;
	ssize_t read;
	size_t rlen;
	char *line;
        int i;

        for (i=0;i<100;i++)
            LOOKUPTABLE[i]=pow(10,-.1*i);

	record rec;

	fp = fopen(argv[1], "r");

	if (fp == NULL)
		exit(EXIT_FAILURE);

	while(1)
	{
		if (getline(&line, &rlen, fp) == -1) break;
		switch(pos%4)
		{
			case 0:
				strcpy(rec.header,line);
				break;
			case 1:
				strcpy(rec.dna,line);
				break;
			case 2:
				strcpy(rec.third,line);
				break;
			case 3:
				strcpy(rec.qua,line);
				break;
		}
		pos++;
		if (!(pos%4))
			recordStats(&rec);
	}

}
