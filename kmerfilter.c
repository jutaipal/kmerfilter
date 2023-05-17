// kmerfilter SCRIPT TO FILTER SHORT READS BASED ON NUMBER OF SHARED KMERS (DEFAULT 8-MER)
// OVERALL DESIGN BY JUSSI TAIPALE, MAY 16 2023
// CODE LARGELY WRITTEN BY CHATGPT4, DEBUGGING AND CORRECTION OF ALGORITHM BY JT
// COMMENTS IN ALL CAPS BY JT, lower case by GPT4
// USAGE: ./kmerfilter file.seq <cutoff>
// GIVING 'n' AS CUTOFF ONLY PRINTS THE NUMBER OF MAX SHARED KMERS BETWEEN THE READ AND PRECEDING PRINTED READS
// GIVING 'h' AS CUTOFF ONLY PRINTS DATA FOR A HISTOGRAM OF SHARED KMER COUNTS
// WILL REMOVE ALL READS THAT CONTAIN MORE THAN cutoff PARTIALLY OR COMPLETELY INDEPENDENT IDENTICAL KMERS
// TO ANY OF THE PREVIOUS READS THAT HAVE BEEN PRINTED (DOES NOT CONSIDER REJECTED READS)
// WHEN A KMER MATCHES ITS REVERSE COMPLEMENT WILL ALWAYS MATCH AS WELL, THIS PAIR IS COUNTED
// AS ONE MATCH. OVERLAPPING KMER MATCHES ARE COUNTED AS SEPARATE EVENTS, AS THEY ARE PARTIALLY INDEPENDENT
// FOR LARGE FILES, ONLY TOP 500,000 READS ARE USED FOR COMPARISON, ALL READS ARE STILL FILTERED
// VERSION 0.11, MAY 17 2023

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// UNCOMMENT BELOW TO FORCE OPENMP ENABLED COMPILATION
#include <omp.h>

// SETS LENGTH OF SUBSEQUENCE (LONGER THAN 8 REQUIRES CHANGE FROM UNSIGNED SHORT TO UNSIGNED LONG INT, ALSO TO SUBROUTINE DNAtoShort)
#define SUBSEQ_LENGTH 8
#define MAX_SAMPLE_SIZE 500000 //LIMITS THE READ POOL THAT EACH READ IS COMPARED TO TO FIRST 500,000 READS TO ALLOW EFFICIENT FILTERING OF LARGE READ FILES

// Global variables for maximum sequence length and number of lines
int MAX_SEQ_LENGTH = 0;
int MAX_LINES = 0;

// Function to calculate maximum sequence length and number of lines
void calculateFileParameters(char *filename) {
    FILE* file = fopen(filename, "r");
    if(file == NULL) {
        printf("Error opening file\n");
        return;
    }

    char line[1000]; // Buffer for reading lines
    while(fgets(line, sizeof(line), file)) {
        MAX_LINES++;
        int len = strlen(line) -1; // -1 to remove newline
        if(len > MAX_SEQ_LENGTH) {
            MAX_SEQ_LENGTH = len;
        }
    }
    fclose(file);
}


// FUNCTION DEFINITIONS (FULL CODE BELOW MAIN PROGRAM)
short int DNAtoShort(char* dna, unsigned short *result);
int countSharedSubsequences(unsigned short *start1, unsigned short *end1, unsigned short *start2, unsigned short *end2);
int getMaxSharedCount(int lineCount, unsigned short **start, unsigned short **end);
int compare(const void* a, const void* b);


// MAIN PROGRAM
int main(int argc, char **argv) {
    
    // CHECKS IF THERE ARE TOO FEW ARGUMENTS
    if(argc != 3) {
        printf("Please input a file name and a shared subsequence limit as arguments. USAGE: ./kmerfilter file.seq <cutoff>\nGiving a cutoff of 'n' does not filter reads, but prints the max number of shared kmers between the read and a preceding read.\nGiving a cutoff of 'h' prints data for a histogram of shared kmer counts. This option can be useful for setting the cut-off.\n");
        return 1;
    }

    calculateFileParameters(argv[1]);
    //printf("\nMax sequence length: %i", MAX_SEQ_LENGTH);
    //printf("\nNumber of lines: %i\n", MAX_LINES);
    
    // IF THERE ARE VERY HIGH NUMBER OF READS, A SAMPLE OF TOP READS IS USED TO FILTER
    if(MAX_LINES > MAX_SAMPLE_SIZE) MAX_LINES = MAX_SAMPLE_SIZE;
    
    // Opens file
    FILE *file = fopen(argv[1], "r");
    
    short int print_max_counts = 0;
    short int print_histogram = 0;
    int sharedSubsequenceLimit = 65536; // INCREASE IF KMER IS LONGER
    if ((argv[2])[0] == 'n') print_max_counts = 1;
    else sharedSubsequenceLimit = atoi(argv[2]);
    if ((argv[2])[0] == 'h') print_histogram = 1;


    if(file == NULL) {
        printf("Error opening file\n");
        return 1;
    }

    // CHECKS IF FILE IS FASTA OR FASTQ
    char first_char = fgetc(file);
    fseek(file, 0, SEEK_SET); // Reset file pointer to the beginning
    int is_fastq = 0;
    if (first_char == '@') is_fastq = 1;
    int is_fasta = 0;
    if (first_char == '>') is_fasta = 1;

    // printf ("fastq: %i, fasta: %i\n", is_fastq, is_fasta);
    short int fastq_read_data_position = 0;
    short int read_status = 0;
    
    // VARIABLES (UNSIGNED SHORT * OK UP TO 8 bp, THEN SHOULD BE CHANGED TO UNSIGNED LONG *)
    unsigned short *results = malloc(sizeof(unsigned short) * 2 * (MAX_LINES+5) * (MAX_SEQ_LENGTH-SUBSEQ_LENGTH+3)+20);
    unsigned short *start[MAX_LINES+2];
    unsigned short *end[MAX_LINES+2];
    char line[MAX_SEQ_LENGTH+2];
    char read_name[MAX_SEQ_LENGTH+2];
    int histogram[MAX_SEQ_LENGTH*2];
    int i;
    for(i = 0; i < MAX_SEQ_LENGTH*2; i++) histogram[i] = 0;
    int lineCount = 0;

    int maxSharedCount = sharedSubsequenceLimit + 1;
    
    
    
    // READ LINE MAIN WHILE LOOP
    while (fgets(line, MAX_SEQ_LENGTH+2, file)) {
        
        line[strcspn(line, "\n")] = '\0';  // remove trailing newline
        // printf("INPUT: %s\n", line);
        // READ NON-SEQUENCE LINES in FASTA and FASTQ
        fastq_read_data_position++;
        if (line[0] == '>' && is_fasta == 1) {
            strncpy (read_name, line, MAX_SEQ_LENGTH); // MOVE THE INPUT LINE TO READ NAME
            continue;
        }
        if (line[0] == '@' && is_fastq == 1) {
            strncpy (read_name, line, MAX_SEQ_LENGTH); // MOVE THE INPUT LINE TO READ NAME
            fastq_read_data_position = 1;
            // printf("READNAME: %s\n", read_name);
            continue;
        }
        if (line[0] == '+' && is_fastq == 1) {
            if(maxSharedCount <= sharedSubsequenceLimit && print_histogram == 0 && read_status == 0) printf("%s\n", line); // PRINTS THE OPTIONAL INFORMATION LINE
            // printf("OPTINAL LINE: %s\n", line);
            continue;
        }
        if (fastq_read_data_position == 4 && is_fastq == 1) {
            if(maxSharedCount <= sharedSubsequenceLimit && print_histogram == 0 && read_status == 0) printf("%s\n", line); // PRINTS THE QUALITY LINE
            fastq_read_data_position = 0;
            // printf("QUALITY: %s\n", optional_information);
            continue;
        }

        
        // SETS START OF THE KMERS FOR THIS READ
        start[lineCount] = &results[2 * lineCount * MAX_SEQ_LENGTH];
        end[lineCount] = start[lineCount];
            
        // GETS KMERS FOR THIS READ
        for(int i = 0, read_status = 0; read_status == 0 && i <= strlen(line) - SUBSEQ_LENGTH; i++) {
            read_status = DNAtoShort(&line[i], end[lineCount]);
            end[lineCount] += 2;
        }
        
        if (read_status == 1) continue; // DISCARDS READS WITH BAD NUCLEOTIDES
          
        //printf("\n");
        //for (i = 0; i < end[lineCount] - start[lineCount]; i++) printf("%i,", *(start[lineCount] + i));
        
        // SORTS THE KMERS
        qsort(start[lineCount], end[lineCount]-start[lineCount], sizeof(unsigned short), compare);
        
        //printf("\n");
        //for (i = 0; i < end[lineCount] - start[lineCount]; i++) printf("%i,", *(start[lineCount] + i));
        //printf("\n");
        
        // GETS MAX SHARED COUNT AGAINST PREVIOUS ALLOWED KMERS
        maxSharedCount = getMaxSharedCount(lineCount, start, end);
        // printf("%i\n", maxSharedCount);
        
        histogram[maxSharedCount]++;  // UPDATES HISTOGRAM
        
        // PRINTS SEQUENCE IF IT DOES NOT SHARE TOO MANY KMERS
        if(maxSharedCount <= sharedSubsequenceLimit) {
            if (is_fastq && print_histogram == 0) printf("%s\n", read_name);
            if (print_max_counts == 0 && print_histogram == 0) printf("%s\n", line);
            if (print_max_counts == 1 && print_histogram == 0) printf("%s\t%i\n", line, maxSharedCount);
            if (lineCount < MAX_LINES) lineCount++; // STOPS INCREMENTING WHEN MAX LINES REACHED
        }
        else {
            // printf("Discarded %s due to %i shared kmers\n", line, maxSharedCount);
        }
    }

    if(print_histogram == 1){
        printf ("Matches\tCount\n");
        for(i = 0; i < MAX_SEQ_LENGTH*2; i++) printf ("%i\t%i\n", i, histogram[i]);
    }
        
    // CLOSES AND EXITS
    fclose(file);
    free(results);

    //
    
    return 0;
}

// SUBROUTINE THAT CONVERTS KMERS TO UNSIGNED SHORT INT (OK UP TO 8-mer, above, change to unsigned long int)
short int DNAtoShort(char* dna, unsigned short *result) {
    result[0] = 0;  // forward
    result[1] = 0;  // reverse
    unsigned short forwardVal = 0;
    unsigned short reverseVal = 0;

    for(int i = 0; i < SUBSEQ_LENGTH; i++) {
        switch(dna[i]) {
            case 'A':
            case 'a':
                forwardVal = 0;
                reverseVal = 3;
                break;
            case 'C':
            case 'c':
                forwardVal = 1;
                reverseVal = 2;
                break;
            case 'G':
            case 'g':
                forwardVal = 2;
                reverseVal = 1;
                break;
            case 'T':
            case 't':
                forwardVal = 3;
                reverseVal = 0;
                break;
            default:
                // printf("Invalid DNA sequence\n");
                return 1; // BAD READ
        }
        result[0] = (result[0] << 2) | forwardVal;
        result[1] = (result[1] >> 2) | (reverseVal << (2 * (SUBSEQ_LENGTH - 1)));
    }
    
    //printf("f%i,r%i\t", result[0], result[1]);
    return 0; // GOOD READ
}

// COMPARE FUNCTION FOR QUICKSORT
int compare(const void* a, const void* b) {
    if (*(unsigned short*)a > *(unsigned short*)b) return 1;
    else if (*(unsigned short*)a < *(unsigned short*)b) return -1;
    else return 0;
}

// NOT CALLED BRUTE FORCE SUBROUTINE THAT COUNTS SHARED SUBSEQUENCES BETWEEN TWO READS (WILL GIVE TWO MATCHES FOR KMER AND ITS REVERSE COMPLEMENT). COUNTS ONLY THE FIRST MATCH. ALLOWS MULTIPLE MATCHES FROM TEST READ KMERS TO PREVIOUS READ KMERS (NOT EASY TO FIX WITHOUT SLOWING DOWN)
int countSharedSubsequences(unsigned short *start1, unsigned short *end1, unsigned short *start2, unsigned short *end2) {
    int count = 0;
    for(unsigned short *ptr1 = start1; ptr1 < end1; ptr1++) {
        for(unsigned short *ptr2 = start2; ptr2 < end2; ptr2++) {
            if(ptr1[0] == ptr2[0]) {
                count++;
                continue; // COUNTS ONLY THE FIRST MATCH
            }
        }
    }
    return count;
}

// SUBROUTINE THAT COUNTS SHARED SUBSEQUENCES BETWEEN TWO READS (WILL GIVE TWO MATCHES FOR KMER AND ITS REVERSE COMPLEMENT). COUNTS ONLY THE FIRST MATCH. REQUIRES SORTED KMERVALUES. MATCHES KMERS 1:1, ANY MATCHING KMERS ARE NOT CONSIDERED AGAIN.
int fastcountSharedSubsequences(unsigned short *start1, unsigned short *end1, unsigned short *start2, unsigned short *end2) {
    int count = 0;
    unsigned short *ptr1 = start1;
    unsigned short *ptr2 = start2;
    
    while(ptr1 < end1 && ptr2 < end2) {
        if(ptr1[0] < ptr2[0]) {
            ptr1++;
        } else if(ptr1[0] > ptr2[0]) {
            ptr2++;
        } else {
            count++;
            ptr1++; // INCREMENTS BOTH POINTERS TO PREVENT
            ptr2++; // MATCHING KMERS FROM BEING CONSIDERED AGAIN
        }
    }
    return count;
}

// SUBROUTINE THAT GETS MAX NUMBER OF SHARED SUBSEQUENCES FOR ALL PRECEDING PASSED READS (OPENMP SUPPORT NOT TESTED)
int getMaxSharedCount(int lineCount, unsigned short **start, unsigned short **end) {
    int maxCount = 0;

    #pragma omp parallel for reduction(max: maxCount)
    for(int i = 0; i < lineCount; i++) {
        int sharedCount = fastcountSharedSubsequences(start[i], end[i], start[lineCount], end[lineCount]);
        if(sharedCount > maxCount) {
            maxCount = sharedCount;
        }
    }
    return maxCount / 2; // DIVIDE BY 2 TO GIVE MAX NUMBER OF PARTIALLY INDEPENDENT MATCHES
}
