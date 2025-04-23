#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/stat.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <wait.h>
#include <sqlite3.h>

/* portul folosit */
#define PORT 2025

/* codul de eroare returnat de anumite apeluri */
extern int errno;


void decode_message_from_client(char* msg, char* decoded_message){

    int len=strlen(msg);
    int cop_i, i;
    int len_decoded_message=0;


    char Base64[65];
	for(i=0;i<26;i++)Base64[i]=i+'A';
	for(i=26;i<52;i++)Base64[i]=i-26+'a';
	for(i=52;i<62;i++)Base64[i]=i-52+'0';
	Base64[62]='+'; Base64[63]='/';Base64[64]='\0';

    int v[4];v[0]=0,v[1]=0,v[2]=0,v[3]=0;int pos;
    for(i = 0; i < len;i++){

        pos=0;
        
        if(msg[i]>='A'&&msg[i]<='Z')pos = msg[i]-'A';else
        if(msg[i]>='a'&&msg[i]<='z')pos = msg[i] - 'a' + 26;else
        if(msg[i]>='0'&&msg[i]<='9')pos = msg[i] - '0' + 52;else
        if(msg[i]=='+')pos=62;else pos=63;

        v[i%4]=pos;
        if(i%4==3){

            cop_i=i;
            int bits[24];for(int j=0;j<24;j++)bits[j]=0;
            for(int j = 0; j < 4; j++){
                int cop=v[j], pos=6*j+5;
                while(cop){// 0 -> 6; 1 -> 12 ; 2-> 18; 3 -> 24 
                    bits[pos]=cop%2;
                    pos--;
                    cop/=2;
                }
            }

            int d1=0,d2=0,d3=0;char c1,c2,c3;
            int putere;
            putere=1;for(int j=7;j>=0;j--){d1+=bits[j]*putere;putere *= 2;}
            putere=1;for(int j=15;j>=8;j--){d2+=bits[j]*putere;putere *= 2;}
            putere=1;for(int j=23;j>=16;j--){d3+=bits[j]*putere;putere *= 2;}

            c1=(char)d1;c2=(char)d2;c3=(char)d3;
            decoded_message[len_decoded_message]=c1; len_decoded_message++;
            decoded_message[len_decoded_message]=c2; len_decoded_message++;
            decoded_message[len_decoded_message]=c3; len_decoded_message++;

            v[0]=0,v[1]=0,v[2]=0,v[3]=0;
        }
    }

        int bits[24];for(int j=0;j<24;j++)bits[j]=0;
            for(int j = 0; j <= len-1-cop_i; j++){
                int cop=v[j], pos=6*j+5;
                while(cop){// 0 -> 6; 1 -> 12 ; 2-> 18; 3 -> 24 
                    bits[pos]=cop%2;
                    pos--;
                    cop/=2;
                }
            }
    int d1=0,d2=0; char c1,c2;
    int putere;
    if(len-1-cop_i==3){
        putere=1;for(int j=7;j>=0;j--){d1+=bits[j]*putere;putere *= 2;}
        putere=1;for(int j=15;j>=8;j--){d2+=bits[j]*putere;putere *= 2;}
        c1=(char)d1;c2=(char)d2;
        decoded_message[len_decoded_message]=c1; len_decoded_message++;
        decoded_message[len_decoded_message]=c2; len_decoded_message++;
    }
    else{
        putere=1;for(int j=7;j>=0;j--){d1+=bits[j]*putere;putere *= 2;}
        c1=(char)d1;
        decoded_message[len_decoded_message]=c1; len_decoded_message++;
    }

    decoded_message[len_decoded_message]='\0';
}

void list(const char *folder, char *output) {
    char command[4096];output[0] = '\0'; 

    snprintf(command, sizeof(command), "ls %s", folder);

    FILE *pipe = popen(command, "r");
    if (pipe == NULL) {
        perror("Eroare la list");
        exit(EXIT_FAILURE);
    }

    char buffer[100];
    
    while (fgets(buffer, sizeof(buffer), pipe) != NULL) {
        strcat(output, buffer);
    }

    if (pclose(pipe) == -1) {
        perror("Eroare la list");
        exit(EXIT_FAILURE);
    }
}

int delfile(char *path){

    if (remove(path) == 0) {
        
        return 0; 
    } else {
        perror("Error deleting file");
        return -1; 
    }
}

void set_length_of_message(int nr_to_transform, char* msg){
    unsigned char* transform=(unsigned char*)&nr_to_transform;
    for(int i=0;i<4;i++)msg[i]=transform[i];msg[4]='\0';
}

int get_length_of_message(int socket){
    char c1,c2,c3,c4;
    read(socket,&c1,1);read(socket,&c2,1);read(socket,&c3,1);read(socket,&c4,1);
    return ((unsigned char)c1 << 24 | (unsigned char)c2 << 16 | (unsigned char)c3 << 8 | (unsigned char)c4);
}

void interpret_message(char *msg, char *msgrasp, char* current_folder, char* copy_folder, char* logged_in, char* USERNAME){//return 1 daca login succes;
                                                                      //return 2 daca signin succes;
                                                                      //return 3 daca exit succes;
                                                                      //return 4 daca logout succes
                                                                      //return 5 daca
    sqlite3 *db;
    sqlite3_stmt *stmt;
    int rc;

    rc = sqlite3_open("accounts.db", &db);
    if (rc){
        strcat(msgrasp, "Eroare: Nu am putut deschide baza de date.");
        return;
    }else printf("S a deschis baza de date.\n");

    if (strncmp(msg, "LOGIN@", 6)==0){

        if(logged_in[0]=='1'){strcpy(msgrasp,"Deja esti logat in cont.");}else{

            char *data = msg+6; 
            char username[50]= {0}, password[4096]={0};
            char *separator=strchr(data, ':');

            int username_length=separator-data;
            strncpy(username, data, username_length);
            username[username_length]= '\0';
            strncpy(password, separator+1, sizeof(password)-1);
            password[sizeof(password)-1] = '\0';

            char *sql = "SELECT COUNT(*) FROM accounts WHERE username = ? AND password = ?";
            if ((rc = sqlite3_prepare_v2(db, sql, -1, &stmt, NULL)) != SQLITE_OK){
                strcat(msgrasp, "Eroare: Nu am putut executa interogarea.");
                sqlite3_close(db);
                return;
            }

            sqlite3_bind_text(stmt,1, username, -1, SQLITE_STATIC);
            sqlite3_bind_text(stmt,2, password, -1, SQLITE_STATIC);

            if ((rc = sqlite3_step(stmt)) == SQLITE_ROW && sqlite3_column_int(stmt, 0)!=0){ 
                strcpy(msgrasp, "Logare realizata cu succes!");
                strcpy(current_folder,"./USERS/");strcat(current_folder,username);
                strcpy(copy_folder,"./COPIES/");strcat(copy_folder,username);
                strcpy(USERNAME,username);
                logged_in[0]='1';
                //strcat(msgrasp,curr_dir);
            }
            else strcpy(msgrasp, "Eroare: Logarea in cont a esuat. Utilizator inexistent / Parola incorecta.");
        
            sqlite3_finalize(stmt);
            
        }

    }else if (strncmp(msg, "SIGNIN@", 7)==0){

        if(logged_in[0]=='1'){strcpy(msgrasp,"Deja esti logat in cont. Delogheaza-te pentru a face un cont nou");return;}else{

            printf("\nse face signin...\n");
            char username[50]={0}, password[4096]={0};//, password_again[50]={0};
            char *data = msg+7;
            char *separator=strchr(data, ':');

            int username_length=separator -data;    //AICEA AI RAMAS, TRB SA SCHIMBI CUM CITESTI MESAJU PRIMIT DE LA CLIENT
            strncpy(username, data, username_length);
            username[username_length]='\0';
            strcpy(password,separator+1);// POATE AICEA MODIFIC CA LA LOGIN

            char *exista_user_sql = "SELECT COUNT(*) FROM accounts WHERE username = ?";
            if ((rc = sqlite3_prepare_v2(db, exista_user_sql, -1, &stmt, NULL))!=SQLITE_OK){
                strcpy(msgrasp, "Eroare: Nu am putut executa verificarea.");
                sqlite3_close(db);
                return;
            }

            sqlite3_bind_text(stmt, 1, username, -1, SQLITE_STATIC);
        
            int user_exists = ((rc = sqlite3_step(stmt))==SQLITE_ROW && sqlite3_column_int(stmt, 0)!=0);
            sqlite3_finalize(stmt);
            if (user_exists) strcpy(msgrasp, "Eroare la creeare cont: utilizatorul deja exista.");
        
            else{
                char *insert_sql = "INSERT INTO accounts (username, password) VALUES (?, ?)";
                if ((rc = sqlite3_prepare_v2(db, insert_sql, -1, &stmt, NULL))!=SQLITE_OK){
                    strcat(msgrasp, "Eroare: Nu am putut executa inserarea.");
                    sqlite3_close(db);
                    return;
                }

                sqlite3_bind_text(stmt, 1, username, -1, SQLITE_STATIC);
                sqlite3_bind_text(stmt, 2, password, -1, SQLITE_STATIC);


                printf("Continutul mesajului pe bucati:\n%s\n%s\n\n",username, password);

    
                if ((rc=sqlite3_step(stmt)) ==SQLITE_DONE){
                
                    //FOLDERU PT USER
                    strcpy(current_folder,"./USERS/");
                    strcat(current_folder,username);//strcat(current_folder,"/");
                    if(mkdir(current_folder,0755)==-1){strcpy(msg,"Nu am putut creea director nou pentru user.");sqlite3_finalize(stmt);return; }

                    //FOLDERU COPIE USER
                    strcpy(current_folder,"./COPIES/");
                    strcat(current_folder,username);//strcat(current_folder,"/");
                    if(mkdir(current_folder,0755)==-1){strcpy(msg,"Nu am putut creea director nou pentru user.");sqlite3_finalize(stmt);return; }

                    strcpy(msgrasp, "Creearea contului a fost realizata cu succes! Sunteti logat in contul nou.");
                    strcpy(current_folder,"./USERS/");strcat(current_folder,username);
                    strcpy(USERNAME,username);
                    logged_in[0]='1';
                }
                else strcpy(msgrasp, "Eroare la creeare cont: problema la insert into.");

                sqlite3_finalize(stmt);
            }
        }
    }else if(strncmp(msg,"exit",4)==0){

        printf("Se inchide conexiunea cu un client.");fflush(stdout);
        sqlite3_close(db);printf("S a inchis baza de date\n");
        exit(0);

    }else if(strncmp(msg,"LOGOUT@",7)==0){

        if(logged_in[0]=='0'){strcpy(msgrasp,"Nu aveti din ce cont sa va delogati.");}else{
            strcpy(msgrasp, "Delogare realizata cu succes.");
            logged_in[0]='0';
            strcpy(current_folder,"./USERS/");
            bzero(USERNAME,sizeof(USERNAME));strcpy(USERNAME,"");
            printf("Un client iese din cont.");fflush(stdout);
        }

    }else if(strncmp(msg,"UPLOAD@",7)==0){//input de forma: UPLOAD@NUME@CONTENT

        if(logged_in[0]=='0'){strcpy(msgrasp,"Nu sunteti logat in cont inca.");}else{

            char copy[4096];
            char nume_fis[100];
            int pos=3;int nr=0;
            for(int i=7;i<strlen(msg);i++)if(msg[i]=='@'){
                pos=i;
                for(int j=7;j<=pos-1;j++){nume_fis[nr]=msg[j];nr++;}
                nume_fis[nr]='\0';
                break;
            }
            strcpy(msgrasp, "Se incarca fisierul in cloud. ");
            strcpy(copy,msg+strlen(nume_fis)+8);

            char cale_fis[4096];
            int fd;
            strcpy(cale_fis,"./USERS/");strcat(cale_fis,USERNAME);strcat(cale_fis,"/");strcat(cale_fis,nume_fis);printf("\n%s\n",cale_fis);
            fd=open(cale_fis,O_WRONLY | O_CREAT | O_TRUNC, 0600);close(fd);
            fd=open(cale_fis,O_RDWR);
            if(write(fd,copy,sizeof(copy))<=0){strcpy(msgrasp,"Eroare la write upload din server(1).");return;};
            close(fd);

            strcpy(cale_fis,"./COPIES/");strcat(cale_fis,USERNAME);strcat(cale_fis,"/");strcat(cale_fis,nume_fis);printf("\n%s\n",cale_fis);
            fd=open(cale_fis,O_WRONLY | O_CREAT | O_TRUNC, 0600);close(fd);
            fd=open(cale_fis,O_RDWR);
            if(write(fd,copy,sizeof(copy))<=0){strcpy(msgrasp,"Eroare la write upload din server(2).");return;};
            close(fd);
            
            strcpy(msgrasp,"Incarcare cu succes!");
        }

    }else if(strncmp(msg,"DOWNLOAD@",9)==0){//raspuns de forma: FILE@nume@content

        if(logged_in[0]=='0'){strcpy(msgrasp,"Nu sunteti logat in cont inca.");}else{
            char nume_fis[100];strcpy(nume_fis,msg+9);//printf("mesaj de la client")
            int pos=3;
            //for(int i=9;i<strlen(msg),i++)if(msg[i]=='@'){pos=i;int nr=0; for(int j=9;j<=i-1;j++)nume_fis[nr]=msg[j],nr++;nume_fis[nr]='\0';break;}
            //strcpy(msgrasp, "Se descarca fisierul din cloud.\n");
            char cale_fis[strlen(current_folder)+100];strcpy(cale_fis,current_folder);strcat(cale_fis,"/");strcat(cale_fis,nume_fis);printf("cale fisier: %s\n",cale_fis);
            int fd=open(cale_fis,O_RDONLY);struct stat sb;fstat(fd,&sb);
            bzero(msgrasp,sizeof(msgrasp));int L=0;
            while((L+=read(fd,msgrasp,sb.st_size))<=0){strcpy(msgrasp,"Eroare la write upload.");};msgrasp[L]='\0';
            close(fd);printf("scos din clould la download: %s\n",msgrasp);
            //strcat(msgrasp,"Descarcare cu succes!");
        }
    }
    else if(strcmp(msg,"LIST@")==0){

        if(logged_in[0]=='0'){strcpy(msgrasp,"Nu sunteti logat in cont inca.");}else{
            char rez[40960];
            list(current_folder,rez);
            strcpy(msgrasp,rez);
        }
    }
    else if(strncmp(msg, "DELFILE@",8)==0){

        if(logged_in[0]=='0'){strcpy(msgrasp,"Nu sunteti logat in cont inca.");}else{
            int rez;
            char cale_fis[100];
            strcpy(cale_fis,current_folder);strcat(cale_fis,"/");
            strcat(cale_fis,msg+8);
            rez=delfile(cale_fis);
            if(rez==0)strcpy(msgrasp,"Fisier sters cu succes!");
            else {strcpy(msgrasp,"Fisier negasit.");}
        }
    }
    else strcpy(msgrasp, "Cerere necunoscuta.");
    sqlite3_close(db);printf("S a inchis baza de date\n");
}

int main ()
{
    struct sockaddr_in server;	// structura folosita de server
    struct sockaddr_in from;
            //mesaj de raspuns pentru client
    int sd;			//descriptorul de socket

    /* crearea unui socket */
    if ((sd = socket (AF_INET, SOCK_STREAM, 0)) == -1)
    {
    	perror ("[server]Eroare la socket().\n");
    	return errno;
    }

    /* pregatirea structurilor de date */
    bzero (&server, sizeof (server));
    bzero (&from, sizeof (from));

    /* umplem structura folosita de server */
    /* stabilirea familiei de socket-uri */
    server.sin_family = AF_INET;
    /* acceptam orice adresa */
    server.sin_addr.s_addr = htonl (INADDR_ANY);
    /* utilizam un port utilizator */
    server.sin_port = htons (PORT);

    /* atasam socketul */
    if (bind (sd, (struct sockaddr *) &server, sizeof (struct sockaddr)) == -1)
    {
    	perror ("[server]Eroare la bind().\n");
    	return errno;
    }

    /* punem serverul sa asculte daca vin clienti sa se conecteze */
    if (listen (sd, 1) == -1)
    {
    	perror ("[server]Eroare la listen().\n");
    	return errno;
    }

    /* servim in mod concurent clientii... */
    while (1)
    {
    	int client;
    	int length = sizeof (from);
    	printf ("[server]Asteptam la portul %d...\n",PORT);
    	fflush (stdout);

    	/* acceptam un client (stare blocanta pina la realizarea conexiunii) */
    	client = accept (sd, (struct sockaddr *) &from, &length);

    	/* eroare la acceptarea conexiunii de la un client */
    	if (client < 0)
    	{
    		perror ("[server]Eroare la accept().\n");
    		continue;
    	}

    	int pid;
    	if ((pid = fork()) == -1) {
    		close(client);
    		continue;
    	} else if (pid > 0) {
    		// parinte
    		close(client);
    		while(waitpid(-1,NULL,WNOHANG));
    		continue;
    	} else if (pid == 0) {

            char msg[4096];
            char msgrasp[4096]="";
            char logged_in[2];
            logged_in[0]='0';logged_in[1]='\0';
            char USERNAME[4096];bzero(USERNAME,sizeof(USERNAME));//USERNAME[0]='\0';
            char current_folder[4096];strcpy(current_folder,"./USERS/");
            char copy_folder[4096];strcpy(copy_folder,"./COPIES/");
            
    		// copil
    		close(sd);
            while(1){
    		    /* s-a realizat conexiunea, se astepta mesajul */
    		    bzero (msg, sizeof(msg));
    		    printf ("[server]Asteptam mesajul...\n");
    		    fflush (stdout);

    		    /* citirea mesajului */

                bzero(msg,sizeof(msg));
                int L;
    		    if ((L=read (client, msg, sizeof(msg))) <= 0)
    		    {
    			    perror ("[server]Eroare la read() de la client.\n");
    			    close (client);	/* inchidem conexiunea cu clientul */
    			    //continue;		/* continuam sa ascultam */
                    exit(0);
    		    }
    		    printf ("[server]Mesajul a fost receptionat...%s,,%ld\n", msg, strlen(msg));fflush(stdout);msg[L-1]='\0';
                //msg[L]='\0';

    		
                //interpret_message(char *msg, char *msgrasp, char* current_folder, char* copy_folder, int& logged_in)



                /******************************************************************************************** */

                bzero(msgrasp,sizeof(msgrasp));
			    interpret_message(msg,msgrasp, current_folder, copy_folder, logged_in, USERNAME);


                //******************************************************************************************* */


    		    /* returnam mesajul clientului */
    		    if (write (client, msgrasp, sizeof(msgrasp)) <= 0)
    		    {
    			    perror ("[server]Eroare la write() catre client.\n");
    			    //continue;		/* continuam sa ascultam */
                    exit(0);
    		    }
    		    else
    			    printf ("[server]Mesajul a fost trasmis cu succes.\n");

            }
    		
    		/* am terminat cu acest client, inchidem conexiunea */
    		close (client);
    		exit(0);
    	}

    }				/* while */
}				/* main */