// cc     fast5ToFq.c   -o fast5ToFq -lhdf5

#include "hdf5.h"
#include <pthread.h>
#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>

pthread_mutex_t DIR_LOCK;
pthread_mutex_t PRINT_LOCK;

char *strstrip(char *s)
{
    size_t size;
    char *end;

    size = strlen(s);

    if (!size)
        return s;

    end = s + size - 1;
    while (end >= s && isspace(*end))
        end--;
    *(end + 1) = '\0';

    while (*s && isspace(*s))
        s++;

    return s;
}

int endsWith(char *str,char *end)
{
    if(strlen(str) >= strlen(end))
    {
        if(!strcmp(str + strlen(str) - strlen(end), end))
        {
            return 1;
        }
    }
    return 0;
}

void readFq (char *path, int D1, int D2)
{
    hid_t       file_id, dataset_id;  /* identifiers */
    long int    ha,hs;

    /* Open an existing file. */
    pthread_mutex_lock(&PRINT_LOCK);
    file_id = H5Fopen(path, H5F_ACC_RDWR, H5P_DEFAULT);

    /* Open an existing dataset. */
    dataset_id = H5Dopen(file_id, "/Analyses/Basecall_1D_000/BaseCalled_template/Fastq", H5P_DEFAULT);
    ha=(long int)H5Dget_offset(dataset_id);
    hs=(long int)H5Dget_storage_size(dataset_id);
    H5Dclose (dataset_id);
    H5Fclose (file_id);
    pthread_mutex_unlock(&PRINT_LOCK);

    int fd=open(path,O_RDONLY);

    char buff[hs];
    lseek(fd,ha,0);
    read(fd,buff,hs);

    close(fd);

    pthread_mutex_lock(&PRINT_LOCK);
    write(1,buff,hs);
    fflush(stdout);
    pthread_mutex_unlock(&PRINT_LOCK);
}

struct params
{
    char *src_path;
    DIR *dp;
};

//void *worker(void *x_void_ptr)
//{
//    struct params* p = (struct params *) x_void_ptr;
//    char buff[16000];
//    struct dirent *ep;
//
//    strcpy(buff,p->src_path);
//    int lenpath = strlen(buff);
//
//    if (buff[lenpath-1]!='/')
//    {
//        buff[lenpath++]='/';
//        buff[lenpath]='\0';
//    }
//
//    while (1)
//    {
//        pthread_mutex_lock(&DIR_LOCK);
//        ep = readdir (p->dp);
//        if (!ep)
//        {
//            pthread_mutex_unlock(&DIR_LOCK);
//            break;
//        }
//        buff[lenpath]='\0';
//        strcat(buff,ep->d_name);
//        pthread_mutex_unlock(&DIR_LOCK);
//
//
//        if (endsWith(buff,".fast5"))
//        {
//            readFq(buff,1,1);
//        }
//    }
//}

int main(int argc, char **argv) 
{
    DIR *dp;
    struct dirent *ep;
    int threads=20;
    struct params p;
    int i=0;

    pthread_t tids[threads];

    if (pthread_mutex_init(&DIR_LOCK, NULL) != 0)
    {
        printf("\n mutex init failed\n");
        return 1;
    }
    if (pthread_mutex_init(&PRINT_LOCK, NULL) != 0)
    {
        printf("\n mutex init failed\n");
        return 1;
    }

    char *line = NULL;
    size_t size;
    while (getline(&line, &size, stdin) != -1)
    {
        strstrip(line);
        readFq(line,1,1);
    }

    //	dp = opendir (argv[1]);
    //
    //	p.dp=dp;
    //	if (dp == NULL)
    //	{
    //		perror ("Couldn't open the directory");
    //		return 1;
    //	}
    //
    //	p.src_path=argv[1];
    //
    //	for (i=0;i<threads;i++)
    //	{
    //		pthread_create(&(tids[i]), NULL, &worker, &p);
    //	}
    //
    //	for (i=0;i<threads;i++)
    //	{
    //		pthread_join(tids[0], NULL);
    //	}
    //
    //
    //	(void) closedir (dp);

    return 0;
}
