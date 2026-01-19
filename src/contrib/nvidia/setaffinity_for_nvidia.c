// Name		: setaffinity_for_nvidia.c
// Authors	: Galen Arnold , Guochun Shi
// Date		: Jan., 2012
// Version 	: 1.1
// Copyright  	: Copyright 2012 
//                The Board of Trustees of the University of Illinois.  
//                All rights reserved.
//
//    compile with -D_GNU_SOURCE to enable CPU_ZERO and similar macros below
#include <stdio.h>
#include <stdlib.h>
#include <sched.h>
#include <unistd.h>
#include <string.h>

int cpulistaffinity(int, int *, int*, char *);
int process_core_string_list(const char* _str, int* list, int* ncores);
int process_core_string_item(const char* str, int* sub_list, int* sub_ncores);

// setaffinity sets the cpu core affinity for the current process by following the
// gpu -> pci_bus -> cpulistaffinity mapping through the /proc and /sys filesystem interfaces
// to the kernel data for the Nvidia driver and the system pci_bus core affinity
//
// This routine should be called early in your program, before large data structures are assigned to
// a (potentially distant) numa node
//
// PPN is mpi processes per node as specified via traditional batch systems, the
// assertion here is that there is a 1:1 mapping of my_rank <-> my_gpu
//
// my_gpu may be passed in, or if it's <0 setaffinity() will choose a gpu using the typical
// rank%PPN logic
//
// for MPI apps, the modulo logic below works with ranks assigned in node-packed order, round-robin
// rank layouts would need different arithmetic to get the correct core <-> gpu mapping
int setaffinity_for_nvidia(int *my_rank, int *PPN, int *my_gpu)
{
		cpu_set_t mycore;
		int mycoreid, myslot;
		int cpu_cores[255];
		int rc_affinity;
		char myhost[80], bus_info[255];
		pid_t mypid;
//fprintf(stdout,"..........................setaffinity_for_nvidia: myrank=%d",*my_rank);

        // intended for the serial case, read MYCUDADEVICE from the environment if present to allow for
        // specifying and utilizing multiple serial gpu jobs per host
        char *mycudadevice;
		mycudadevice= getenv("MYCUDADEVICE");
		// if set, MYCUDADEVICE will override the passed value of my_gpu 
		if (mycudadevice != NULL)
				*my_gpu= atoi(mycudadevice);

	    // The MPI_PPN environment setting is targeted at MPI applications
        // read MPI_PPN processes-per-node from the environment so that the code doesn't
        // need to be rebuilt for changes to that common batch environment parameter
        char * myPPN;
		myPPN= getenv("MPI_PPN");
		// if MPI_PPN is set, it will override the passed value of PPN
		if (myPPN != NULL)
				*PPN= atoi(myPPN);

		gethostname(myhost, 80);
		mypid= getpid();
		if ((*my_gpu) < 0 ) *my_gpu= (*my_rank) % (*PPN);
//fprintf(stdout,"..........................setaffinity_for_nvidia: myrank=%d; myPPN=%s; PPN=%d; my_gpu=%d; myost=%s  \n",
//                              *my_rank, myPPN, *PPN, *my_gpu, myhost );

		int ncores =255;

		// call cpulistaffinity to track down the cpu_cores for the provided my_gpu by querying
		// the appropriate files under /proc/driver/nvidia/ and /sys/class/pci_bus/
		cpulistaffinity(*my_gpu, cpu_cores, &ncores, bus_info);
		// CPU_ZERO and similar are macros , see: man CPU_ZERO for more information
		CPU_ZERO(&mycore);
		myslot= (*my_gpu) % ncores ;   // ncores is num elements in cpu_cores[]
		CPU_ZERO(&mycore);
		CPU_SET( cpu_cores[myslot], &mycore );

		// set and check the core affinity, report back the current layout
		rc_affinity= sched_setaffinity(mypid, sizeof(cpu_set_t), &mycore);
		if (rc_affinity != 0) perror("sched_setaffinity");
		rc_affinity= sched_getaffinity(mypid, sizeof(cpu_set_t), &mycore);
		if (rc_affinity != 0) perror("sched_getaffinity");
		if (CPU_ISSET(cpu_cores[myslot], &mycore ))
			mycoreid= cpu_cores[myslot];
	    // debug
	    //	fprintf(stdout, "host %s mpi_rank %d gpu %d cpu_core %d ncores %d myslot %d\n",
	    //			myhost, my_rank, my_gpu,  cpu_cores[myslot], ncores, myslot);
printf("setaffinity_for_nvidia: host %s mpi_rank %d process %d gpu %d cpu_core %d pci_bus %s\n",
                                myhost, *my_rank, mypid, *my_gpu, mycoreid, bus_info);

		// the diagnostic here can be redirected to stderr if desired by swapping it for stdout
		// or comment the fprintf() lines if not desired
		// for validation, compare with nvidia-smi , "numactl --hardware", and
		// /sys/class/pci_bus/*/cpulistaffinity
		fprintf(stdout, "host %s mpi_rank %d process %d gpu %d cpu_core %d pci_bus %s\n",
				myhost, *my_rank, mypid, *my_gpu, mycoreid, bus_info);
		return(*my_gpu);
}

// cpulistaffinity() makes the association between the numbered gpu device "my_gpu" and
// the cpu cores associated with it by following the path name of the pci_bus listed
// under the nvidia driver in /proc .
int cpulistaffinity(int my_gpu, int *cpu_cores, int* ncores,  char * bus_info)
{
	FILE *nvidia_info, *pci_bus_info;
	size_t nbytes = 255;
	//int core3, core4; // fillers for sscanf()
	char *my_line;
	char nvidia_info_path[255], pci_bus_info_path[255];

	// the nvidia driver populates this path for each gpu
	sprintf(nvidia_info_path,"/proc/driver/nvidia/gpus/%d/information", my_gpu);
	nvidia_info= fopen(nvidia_info_path,"r");
	if (nvidia_info == NULL)
	{
		perror(nvidia_info_path);
		exit(-1);
	}

	my_line= (char *) malloc(nbytes +1);
	if (my_line == NULL)
	{ perror("error allocating memory for my_line string"); exit(-1); }

	while (!feof(nvidia_info)) // reading lines in the information file for the nvidia driver
	{
			if ( -1 == getline(&my_line, &nbytes, nvidia_info))
				break;
			else
			{ // the first 7 char of the Bus Location will lead to the corresponding
			  // path under /sys/class/pci_bus/  , cpulistaffinity showing cores on that
			  // bus is located there
				if ( 1 == sscanf(my_line,"Bus Location: %s", bus_info ))
				{
					sprintf(pci_bus_info_path,"/sys/class/pci_bus/%.7s/cpulistaffinity",
							bus_info);
				}
			}
	}
	// open the cpulistaffinity file on the pci_bus for "my_gpu"
	pci_bus_info= fopen(pci_bus_info_path,"r");
	if (pci_bus_info == NULL)
	{
		perror(pci_bus_info_path);
		exit(-1);
	}
	while (!feof(pci_bus_info))
	{
		if ( -1 == getline(&my_line, &nbytes, pci_bus_info))
			break;
		else
		{
			int rc = process_core_string_list(my_line, cpu_cores, ncores);
                        if(rc < 0){
                          printf("ERROR: processing the line (%s) failed\n", my_line);
                          return  -1;
                        }
// debug, devel output
//			int i;
//			printf("gpu %d # cpu cores: %d", my_gpu, *ncores);
//			printf(" list :\t");
//			for( i=0;i < *ncores; i++){
//			  printf("%d\t", cpu_cores[i]);
//			}
//			printf("\n");
		}
	}
	if (my_line) free(my_line);
	return(0);
}


int process_core_string_list(const char* _str, int* list, int* ncores)
{
  /* The input string @str should be separated by comma, and each item can be 
   * either a number or a range (see the comments in process_core_string_item 
   * function)
   *
   */

  if(_str == NULL || list == NULL || ncores == NULL
     || *ncores <= 0){
    printf("ERROR: Invalid arguments in function %s\n", __FUNCTION__ );
    return  -1;
  }

  char str[256];
  strncpy(str, _str, sizeof(str));

  int left_space = *ncores;
  int tot_cores = 0;

  char* item = strtok(str, ",");
  if(item == NULL){
    printf("ERROR: Invalid string format(%s)\n", str);
    return -1;
  }
 
  do {
    int sub_ncores = left_space;
    int* sub_list = list + tot_cores;

    int rc = process_core_string_item(item, sub_list, &sub_ncores);
    if(rc <0){
      printf("ERROR: processing item(%s) failed\n", item);
      return -1;
    }

    tot_cores += sub_ncores;
    left_space -= sub_ncores;

    item = strtok(NULL, ",");
  }while( item != NULL);

  *ncores = tot_cores;
  return 0;
}


int process_core_string_item(const char* str, int* sub_list, int* sub_ncores)
{
  /* assume the input format is one of the following two
   * 1. a number only, e.g. 5
   * 2. a range, e.g 4-6, which means three numbers 4,5,6
   * return a list of numbers in @sub_list and and the total numbers
   * in @sub_ncores
   */
  int i;
  if(str == NULL || sub_list == NULL || sub_ncores == NULL ||
     *sub_ncores <= 0){
    printf("ERROR: Wrong parameters in function %s!\n", __FUNCTION__);
    return -1;
  }

  if(strstr(str, "-") != NULL){
    //a range
    int low_core, high_core;
    if (sscanf(str,"%d-%d",&low_core, &high_core) != 2){
      printf("ERROR: range scan failed\n");
      return -1;
    }
    if(*sub_ncores <  high_core-low_core +1){
      printf("ERROR: not enough space in sub_list\n");
      return -1;
    }

    for(i = 0; i < high_core-low_core +1; i++){
      sub_list[i] = i + low_core;
    }
    *sub_ncores =  high_core - low_core +1;

  }else{
    //a number
    int core;
    if (sscanf(str, "%d", &core) != 1){
      printf("ERROR: wrong format for core number\n");
      return -1;
    }
    sub_list[0] = core;
    *sub_ncores   =1;
  }
  return 0;
}
