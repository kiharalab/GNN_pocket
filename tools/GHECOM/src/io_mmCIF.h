/*
 <io_mmCIF.h>


*/

#define HASH_SIZE 11
#define MAX_LENGTH_LINE 1000
#define MAX_LENGTH_WORD 50000




struct KEY_INT{
  /* "struct KEY_INT " is mainly used for hash_table for connecting 'column_name' and 'column_number'. */
  char *key;             /* (malloc later) */
  int  int_value;
  struct KEY_INT *next;
};


struct ROWLINE{
  /* "struct ROWLINE" is for ROWLINE for the "TABLE". */

  char **data; /* data[Ncolumn][length_of_value_string] */
  struct ROWLINE *next;
};



struct WORD{
  /* for  a list of strings (words) */  
  char *word;     /* (malloc later) */
  struct WORD *next;
};



struct TABLE{
/*
             | column0 | column1 | ........ | column[Ncolumn-1] |
--------------------------------------------------------------
0    rowline |         |         | .........|                   |   
1    rowline |         |         | .........|                   |   
:
:
Nrow rowline |         |         | .........|                   |   
--------------------------------------------------------------
*/

  char *table_name;   /* ID string for the TABLE (malloc later) */

  int  Ncolumn;       /* number of columns */
  int  Nrow;          /* number of rows    */
 
  struct WORD head_column_name;               /* head of list of column name */
  struct KEY_INT  colnum_hash_tab[HASH_SIZE]; /* hash_table for connecting 'column_name' and 'column_number'*/

  struct ROWLINE head_rowline; /* head of linked list of ROWLINE */
  
  struct WORD  temp_head_word; /* temporary WORD linked list only for "ONE-ROW" table. */ 

  struct TABLE *next; /* used for the linked list starting from the BLOCK */
};


struct KEY_TABLE{
  char   *key;    /* (malloc later) */
  struct TABLE *table_value;
  struct KEY_TABLE *next;
};


struct BLOCK{
  /* The "BLOCK" contains the all of the data in the mmCIF file. */
  /* It contains the linked list of several TABLEs. */
  char   *block_name; /* ID for the BLOCK (malloc later) */ 
  int    Ntable;
  struct TABLE head_table; /* head of TABLE linked list  */
  struct KEY_TABLE table_hash_tab[HASH_SIZE];
};


/* Functions (GLOBAL) */
extern void read_mmCIF_file_into_BLOCK();
extern void read_mmCIF_file_into_WORDs();
extern void translate_mmCIF_WORDs_into_BLOCK_data();
extern int data_malloc_from_rowline();
extern char* data_from_rowline();
extern int  number_from_column();
extern int  Integer_from_key();
extern struct TABLE *TABLE_from_key();
extern int free_WORDs();
extern void show_WORDs();
extern void free_BLOCK();
extern void split_into_two_words();

