#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <arpa/inet.h>
#include <errno.h>
#include <stdbool.h>
#include <openssl/sha.h>
#include <fcntl.h>
    

char* select_file() {
    
    const char* command="zenity --file-selection --title=\"Select a file to upload\"";
    char path[1024];
    FILE* fp;
    
    fp=popen(command, "r");
    if (fp==NULL){perror("Eroare la zenity");return NULL;}

    if(fgets(path,sizeof(path),fp) == NULL){pclose(fp);return NULL;}
    pclose(fp);

    size_t len=strlen(path);
    if(len >0 &&path[len-1]=='\n') path[len -1]='\0';

    char* result=(char*)malloc(len);
    if(result != NULL) strcpy(result, path);
    
    return result;
}
   

void encode_file_content_to_send(const char *filename, char *output) {

    char base64_table[65];
	for(int i=0;i<26;i++)base64_table[i]=i+'A';
	for(int i=26;i<52;i++)base64_table[i]=i-26+'a';
	for(int i=52;i<62;i++)base64_table[i]=i-52+'0';
	base64_table[62]='+'; base64_table[63]='/';base64_table[64]='\0';

    FILE *file = fopen(filename, "rb"); 
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);

    unsigned char *buffer = (unsigned char *)malloc(file_size);
    if (buffer == NULL) {
        perror("Memory allocation failed");
        fclose(file);
        return ;
    }

    fread(buffer, 1, file_size, file);
    fclose(file);

    int output_len = 0;
    int i;
    for (i = 0; i < file_size; i += 3) {
        unsigned char byte0 = buffer[i];
        unsigned char byte1; if(i + 1 < file_size) byte1= buffer[i + 1];else byte1= 0;
        unsigned char byte2; if(i + 2 < file_size) byte2= buffer[i + 2];else byte2= 0;

        output[output_len++] = base64_table[byte0 >> 2];
        output[output_len++] = base64_table[((byte0 & 0x03) << 4) | (byte1 >> 4)];
        if(i + 1 < file_size)output[output_len++] = base64_table[((byte1 & 0x0F) << 2) | (byte2 >> 6)]; else output[output_len++] =  '=';
        if(i + 2 < file_size)output[output_len++] = base64_table[byte2 & 0x3F]; else output[output_len++]='=';
    }
    output[output_len] = '\0';
    free(buffer);
}

void decode_message_from_server(char* msg, char* decoded_message){

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

void set_length_of_message(int nr_to_transform, char* msg){
    unsigned char* transform=(unsigned char*)&nr_to_transform;
    for(int i=0;i<4;i++)msg[i]=transform[i];msg[4]='\0';
}

int get_length_of_message(int socket){
    char c1,c2,c3,c4;
    read(socket,&c1,1);read(socket,&c2,1);read(socket,&c3,1);read(socket,&c4,1);
    return ((unsigned char)c1 << 24 | (unsigned char)c2 << 16 | (unsigned char)c3 << 8 | (unsigned char)c4);
}

void interpret_response(char* rasp, char* command){


    printf("%s\n",rasp);

    if(strncmp(command, "download ",9)==0){

        char decoded_text[4096];
        decode_message_from_server(rasp,decoded_text);

        printf("\ndecoded text: %s\n",decoded_text);fflush(stdout);

        int fd=open(command+9,O_WRONLY | O_CREAT | O_TRUNC, 0600);
        write(fd,decoded_text,strlen(decoded_text));
        close(fd);
        printf("Fisierul a fost descarcat cu succes!\n");
    }
}

int main(int argc, char *argv[]) {

    //SA NU UIT DE REDUNDANTA SI ENCRIPTIE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    int sd;                        // descriptorul de socket
    struct sockaddr_in server;     // structura folosita pentru conectare
    char msg[4096];bzero(msg,sizeof(msg));                 // mesajul trimis
    char command[4096];             // comanda utilizatorului
    int port;
    int nrbytes = 0;

    /* verificare argumente */
    if (argc != 3) {
        printf("[client] Sintaxa: %s <adresa_server> <port>\n", argv[0]);
        return -1;
    }

    /* stabilim portul */
    port = atoi(argv[2]);

    /* cream socketul */
    if ((sd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
        perror("[client] Eroare la socket().\n");
        return errno;
    }

    /* configuram structura pentru server */
    server.sin_family = AF_INET;
    server.sin_addr.s_addr = inet_addr(argv[1]);
    server.sin_port = htons(port);

    /* conectare la server */
    if (connect(sd, (struct sockaddr *)&server, sizeof(server)) == -1) {
        perror("[client] Eroare la connect().\n");
        return errno;
    }

    printf("Bine ati venit pe GabiCloud!\n");
    printf("Pentru logare in cont, scrieti 'login'\n");
    printf("Pentru crearea unui cont nou, scrieti 'cont nou':\n");
    printf("Pentru iesire din aplicatie, scrieti 'exit'\n");
    printf("Pentru iesire din cont, scrieti 'logout'\n");
    printf("Pentru incarcare fisier, scrieti 'upload'\n");
    printf("Pentru descarcare fisier, scrieti 'download nume_fisier'\n");

    printf("Pentru vizualizare fisiere / foldere accesibile din directorul curent, scrieti 'list'\n");
    //printf("Pentru creeare folder nou in directorul curent, scrieti 'folder nume_folder'\n");
    //printf("Pentru schimbare director curent, scrieti 'cd nume_director'\n");
    printf("Pentru stergere fisier din directorul curent, scrieti 'delfile nume_fisier'");
    //printf("Pentru stergere un folder accesibil din directorul curent, scrieti 'delfolder nume_folder'\n");

    /* Citirea comenzii de la utilizator */
    
    while (true) {

        bzero(msg, sizeof(msg));
        while(true) {
            printf("\nScrieti comanda: "); fflush(stdout);
            bzero(command, sizeof(command));
            nrbytes=read(0,command,sizeof(command));
            command[nrbytes-1] = '\0';
            //printf("comanda este: %s\n",command);

            if(strcmp(command, "login") == 0 || strcmp(command, "cont nou") == 0 || strcmp(command, "exit") == 0 || strcmp(command, "logout") == 0 || strcmp(command, "upload") == 0 || strncmp(command, "download ",9) == 0 || strcmp(command, "list") == 0 || strncmp(command, "folder ",7) == 0 || strncmp(command, "cd ",3) == 0 || strncmp(command, "delfile ",8) == 0 || strncmp(command, "delfolder ",10) == 0) break;
            
            printf("Comanda nerecunoscuta. Scrieti doar 'login' sau 'cont nou' sau 'exit'.\n");fflush(stdout);
        }

        if(strcmp(command, "login") == 0){

            strcat(msg, "LOGIN@");
            char username[50], password[2*SHA256_DIGEST_LENGTH+1];

            printf("User: ");fflush(stdout);
            nrbytes = read(0,username, sizeof(username));
            username[nrbytes-1]='\0';
            strcat(msg, username);strcat(msg, ":");
            printf("Parola: ");fflush(stdout);
            nrbytes = read(0,password,sizeof(password));
            password[nrbytes-1]='\0';

            unsigned char hash[SHA256_DIGEST_LENGTH];
            SHA256((unsigned char*)password, strlen(password),hash);
            strcpy(password,"");
            for(int i =0;i < SHA256_DIGEST_LENGTH;i++)sprintf(password+(2*i),"%02x",hash[i]);
            password[SHA256_DIGEST_LENGTH*2]='\0';/*password[0]='@';*/
            strcat(msg, password);

        } else if(strcmp(command, "cont nou") == 0){

            strcat(msg, "SIGNIN@");
            char username[50], password[2*SHA256_DIGEST_LENGTH+1], password_again[2*SHA256_DIGEST_LENGTH+1];

            printf("User: ");fflush(stdout);
            nrbytes = read(0,username, sizeof(username));
            username[nrbytes-1]= '\0';
            strcat(msg, username);strcat(msg, ":");

            while(true){
                printf("Parola: ");fflush(stdout);
                nrbytes = read(0,password,sizeof(password));
                password[nrbytes-1]='\0';

                printf("Repetati parola: "); fflush(stdout);
                nrbytes = read(0,password_again,sizeof(password_again));
                password_again[nrbytes-1]='\0';

                if (strcmp(password, password_again)==0) {
                    unsigned char hash[SHA256_DIGEST_LENGTH];
                    SHA256((unsigned char*)password, strlen(password),hash);
                    strcpy(password,"");
                    for(int i =0;i < SHA256_DIGEST_LENGTH;i++)sprintf(password+(2*i),"%02x",hash[i]);
                    password[SHA256_DIGEST_LENGTH*2]='\0';/*password[0]='@';*/strcat(msg, password);break;
                }
                else printf("Parolele nu sunt identice! Incercati din nou.\n");
            }
            
        }else if(strcmp(command, "exit")==0){

            printf("Se inchide conexiunea.\n");fflush(stdout);
            strcat(msg, "exit");
            if(write(sd, msg, sizeof(msg))<=0) {
                perror("Eroare la write spre server.\n");
                close(sd);
            }
            exit(0);

        }else if(strcmp(command, "logout") == 0){

            printf("Iesim din cont.\n");fflush(stdout);
            strcat(msg, "LOGOUT@");

        }else if(strcmp(command, "upload") == 0){

            char* file_to_send=select_file();
            printf("Se trimite fisierul: %s\n",file_to_send);
            printf("Se incarca fisier / folder\n");fflush(stdout);

            struct stat sb;stat(file_to_send,&sb);

            strcat(msg, "UPLOAD@");
            for(int i=strlen(file_to_send);i>=0;i--)if(file_to_send[i]=='/'){strcat(msg,file_to_send+i+1);break;}strcat(msg,"@");
            char* buffer=(char*)malloc(4*sb.st_size+10);encode_file_content_to_send(file_to_send,buffer);strcat(msg,buffer);
            //int fd=open(file_to_send,O_RDONLY);read(fd,buffer,sizeof(buffer));buffer[sb.st_size-1]='\0';
            free(buffer);//0...n-1 ->n ch ||||  sb.st_size=n+1
            

        }else if(strncmp(command, "download ",9) == 0){

            printf("Se descarca fisier\n%s\n",command);fflush(stdout);
            strcpy(msg, "DOWNLOAD@");
            strcat(msg,command+9);
        }
        else if(strcmp(command,"list")==0){

            strcat(msg,"LIST@");
        }
        else if(strncmp(command, "folder ",7) == 0){

            strcat(msg,"FOLDER@");
            strcat(msg,command+7);

        }
        else if(command[0]=='c' && command[1]=='d'){
           
            strcat(msg,"CD@");
            strcat(msg,command+3);         
        }
        else if(strncmp(command, "delfile ", 8)==0){
            
            strcat(msg,"DELFILE@");
            strcat(msg,command+8);
        }
        else if(strncmp(command, "delfolder ", 10)==0){

            
            strcat(msg,"DELFOLDER@");
            strcat(msg,command+10);
        }

        /* trimiterea mesajului la server */
        printf("MESAJUL TRIMIS SPRE SERVER ESTE: %s\n", msg);



        if(write(sd, msg, sizeof(msg)) <= 0) {
            perror("[client] Eroare la write() spre server.\n");
            close(sd);
            continue;
            //return errno;
        }

        /* citirea raspunsului de la server */
        bzero(msg, sizeof(msg));
        if ((read(sd, msg, sizeof(msg))) <= 0) {
            perror("[client] Eroare la read() de la server.\n");
            close(sd);
            continue;
            //return errno;
        }
        //msgrasp[L-1]='\0';

        /* gestionarea raspunsului */
        interpret_response(msg, command);
    }

    /* inchidem conexiunea */
    close(sd);
    return 0;
}
