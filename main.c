#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
    unsigned char B,G,R;

}pixel;

 typedef struct
 {
     unsigned char padding;
     unsigned int size, width, height;
     unsigned int *header;

 }BMP_data;

typedef struct
{
    int digit;
    unsigned short elim;
    int line, col;
    unsigned long size;
    double corr;///corelatia

}window_list;

 int check_memory(const void *data)
 {
     if(data==NULL)
     {
         printf("ERROR: MEMORY ALLOCATION");
         return 1;
     }
     return 0;
 }
BMP_data * get_data(FILE *in)
{
    BMP_data *data;
    data=(BMP_data *)malloc(sizeof(BMP_data));
    if(check_memory(data)) return NULL;
    data->header=(unsigned int *)malloc(54*sizeof(unsigned char));
    if(check_memory(data->header))return NULL;
    fread(data->header,54,1,in);
    fseek(in,2,SEEK_SET);
    fread(&(data->size),sizeof(unsigned int),1,in);
    fseek(in,18,SEEK_SET);
    fread(&(data->width),sizeof(unsigned int),1,in);
    fread(&(data->height),sizeof(unsigned int),1,in);
    if(data->width % 4 != 0)
        data->padding = 4 - (3 * data->width) % 4;
    else
        data->padding = 0;
    return data;
}
pixel *load_in(FILE *in)
{
    BMP_data *info;
    info=(BMP_data*)malloc(sizeof(BMP_data));
    if(check_memory(info))return NULL;
    info=get_data(in);

    pixel *image;
    image=(pixel *)malloc(sizeof(pixel)*info->height*info->width);
    if(check_memory(image))return NULL;

    fseek(in,54,SEEK_SET);
    int i,j,k=0;
    for(i=info->height-1;i>=0;i--)
    {
        for(j=0;j<info->width;j++)
        {
            fread(&image[i * info->width + j].B, 1, 1, in);
            fread(&image[i * info->width + j].G, 1, 1, in);
            fread(&image[i * info->width + j].R, 1, 1, in);
        }
        fseek(in,info->padding,SEEK_CUR);
    }

    return image;
}
void load_out(FILE *out, pixel *linear, BMP_data *info)
{
    fwrite(info->header,54,1,out);
    int i,k=1,j;
    int a=0;
    for(i=info->height-1;i>=0;i--)
    {
        for(j=0;j<info->width;j++)
        {
            fwrite(&linear[i * info->width + j].B,1,1,out);
            fwrite(&linear[i * info->width + j].G,1,1,out);
            fwrite(&linear[i * info->width + j].R,1,1,out);
        }
        fwrite(&a,sizeof(int),info->padding,out);
    }
}
unsigned int * xorshift32(unsigned int seed,unsigned int width,unsigned int height)
{
    unsigned int *random;
    random=(unsigned int*)malloc(2*width*height*sizeof(unsigned int));
    if(check_memory(random))return NULL;
    unsigned int r, i;
    r=seed;
    random[0]=r;
    for(i=1;i<=2*width*height-1;i++)
    {
        r=r^r<<13;
        r=r^r>>17;
        r=r^r<<5;
        random[i]=r;
    }
    return random;
}
void durstenfeld(unsigned int size, unsigned int *sequence, unsigned int *random)
{
    unsigned int i,aux,j,k=1;
    for(i=0;i<size;i++)
        sequence[i]=i;
    for(i=size-1;i>0;i--)
    {
        j=random[k]%(i+1);
        aux=sequence[j];
        sequence[j]=sequence[i];
        sequence[i]=aux;
        k++;
    }
}
pixel * apply_perm(pixel *image, unsigned int *permutation,unsigned int size)
{
    pixel *new_image=(pixel*)malloc(sizeof(pixel)*size);
    if(check_memory(new_image))return NULL;
    unsigned int i;
    for(i=0;i<size;i++)
    {
        new_image[permutation[i]]=image[i];
    }
    return new_image;
}
pixel px_xor_nr(pixel px,unsigned int n)
{
    pixel new_px;
    new_px.B=px.B^(n&255);
    new_px.G=px.G^((n>>8)&255);
    new_px.R=px.R^((n>>16)&255);
    return new_px;
}
pixel px_xor_px(pixel px1,pixel px2)
{
    pixel new_px;
    new_px.B=px1.B^px2.B;
    new_px.G=px1.G^px2.G;
    new_px.R=px1.R^px2.R;
    return new_px;
}
pixel *ciphered_image(pixel *image, unsigned int SV, unsigned int *random,unsigned int size)
{
    pixel *new_image=(pixel *)malloc(sizeof(pixel)*size);
    if(check_memory(new_image))return NULL;
    int i;
    new_image[0]=px_xor_nr(px_xor_nr(image[0],SV),random[size]);
    for(i=size+1;i<=2*size-1;i++)
    {
        new_image[i-size]=px_xor_nr(px_xor_px(new_image[i-size-1],image[i-size]),random[i]);
    }
    return new_image;
}
void encrypt(const char *initial_bmp,const char *final_bmp,const char *key)
{
    FILE *ibmp=fopen(initial_bmp,"rb");
    FILE *fbmp=fopen(final_bmp,"wb");
    FILE *sk=fopen(key,"r");
    if(ibmp==NULL||fbmp==NULL||sk==NULL)
    {
        printf("ERROR:FILE");
        return;
    }
    unsigned int seed,sv;
    fscanf(sk,"%u",&seed);
    fscanf(sk,"%u",&sv);
    ///0.Liniarizare
    BMP_data *info;
    info=(BMP_data*)malloc(sizeof(BMP_data));
    if(check_memory(info))return;
    info=get_data(ibmp);
    pixel *initial_image;
    initial_image=(pixel*)malloc(sizeof(pixel)*info->height*info->width);
    if(check_memory(initial_image))return;
    initial_image=load_in(ibmp);
    ///1.Generare cu XORSHIFT
    unsigned int *random;
    random=(unsigned int*)malloc(2*info->height*info->width*sizeof(unsigned int));
    if(check_memory(random))return;
    random=xorshift32(seed,info->width,info->height);
    ///2.Generare permutare
    unsigned int *permutation;
    permutation=(unsigned int*)malloc(info->height*info->width*sizeof(unsigned int));
    if(check_memory(permutation))return;
    durstenfeld(info->height*info->width,permutation,random);
    ///3.Aplicare permutare
    pixel *interm_image=(pixel*)malloc(sizeof(pixel)*info->height*info->width);
    if(check_memory(interm_image))return;
    interm_image=apply_perm(initial_image,permutation,info->height*info->width);
    ///4.Criptare finala
    pixel *final_image=(pixel*)malloc(sizeof(pixel)*info->height*info->width);
    if(check_memory(final_image))return;
    final_image=ciphered_image(interm_image,sv,random,info->height*info->width);
    ///5.Salvare
    load_out(fbmp,final_image,info);

    fclose(ibmp); fclose(fbmp); fclose(sk);
    free(info); free(initial_image);free(random);
    free(permutation); free(interm_image); free(final_image);

}
unsigned int *inverse_perm(unsigned int *permutation, unsigned int size)
{
    unsigned int *inversed=(unsigned int*)malloc(size*sizeof(unsigned int));
    if(check_memory(inversed))return NULL;
    unsigned int i;
    for(i=0;i<size;i++)
    {
        inversed[permutation[i]]=i;
    }
    return inversed;
}
pixel *deciphered_image(pixel *image, unsigned int SV, unsigned int *random,unsigned int size)
{
    pixel *new_image=(pixel *)malloc(sizeof(pixel)*size);
    if(check_memory(new_image))return NULL;
    unsigned int i;
    new_image[0]=px_xor_nr(px_xor_nr(image[0],SV),random[size]);
    for(i=size+1;i<=2*size-1;i++)
    {
        new_image[i-size]=px_xor_nr(px_xor_px(image[i-size-1],image[i-size]),random[i]);
    }
    return new_image;
}
void decrypt(const char *initial_bmp,const char *final_bmp,const char *key)
{
    FILE *ibmp=fopen(initial_bmp,"rb");
    FILE *fbmp=fopen(final_bmp,"wb");
    FILE *sk=fopen(key,"r");
    if(ibmp==NULL||fbmp==NULL||sk==NULL)
    {
        printf("ERROR:FILE");
        return;
    }
    unsigned int seed,sv;
    fscanf(sk,"%u",&seed);
    fscanf(sk,"%u",&sv);
    ///0.Liniarizare
    BMP_data *info;
    info=(BMP_data*)malloc(sizeof(BMP_data));
    if(check_memory(info))return;
    info=get_data(ibmp);
    pixel *initial_image;
    initial_image=(pixel*)malloc(sizeof(pixel)*info->height*info->width);
    if(check_memory(initial_image))return;
    initial_image=load_in(ibmp);
    ///1.Generare cu XORSHIFT
    unsigned int *random;
    random=(unsigned int*)malloc(2*info->height*info->width*sizeof(unsigned int));
    if(check_memory(random))return;
    random=xorshift32(seed,info->width,info->height);
    ///2.Generare permutare si inversarea ei
    unsigned int *permutation;
    permutation=(unsigned int*)malloc(info->height*info->width*sizeof(unsigned int));
    if(check_memory(permutation))return;
    durstenfeld(info->height*info->width,permutation,random);
    unsigned int *inversed_permutation;
    inversed_permutation=inverse_perm(permutation,info->width*info->height);
    ///3.Inversarea substitutiei initiale
    pixel *interm_image=(pixel*)malloc(sizeof(pixel)*info->height*info->width);
    if(check_memory(interm_image))return;
    interm_image=deciphered_image(initial_image,sv,random,info->height*info->width);
    ///4.Permutarea finala a pixelilor
    pixel *final_image=(pixel*)malloc(sizeof(pixel)*info->height*info->width);
    if(check_memory(final_image))return;
    final_image=apply_perm(interm_image,inversed_permutation,info->height*info->width);
    ///5.Salvare
    load_out(fbmp,final_image,info);

    fclose(ibmp); fclose(fbmp); fclose(sk);
    free(info); free(initial_image);free(random);
    free(permutation); free(interm_image); free(final_image);
}
void chi_test(const char *name, const char *type)
{
    FILE *in=fopen(name,"rb");
    if(in==NULL)
    {
        printf("ERROR:FILE");
        return;
    }
    BMP_data *info;
    info=(BMP_data*)malloc(sizeof(BMP_data));
    if(check_memory(info))return;
    info=get_data(in);
    pixel *image;
    image=(pixel*)malloc(sizeof(pixel)*info->height*info->width);
    if(check_memory(image))return;
    image=load_in(in);

    unsigned int *R=(unsigned int*)calloc(256,sizeof(unsigned int));
    unsigned int *G=(unsigned int*)calloc(256,sizeof(unsigned int));
    unsigned int *B=(unsigned int*)calloc(256,sizeof(unsigned int));
    if(check_memory(R)||check_memory(G)||check_memory(B))return;
    unsigned int i;
    for(i=0;i<info->width*info->height;i++)
    {
        R[image[i].R]++;
        G[image[i].G]++;
        B[image[i].B]++;
    }
    double *chi_result=(double*)calloc(3,sizeof(double));
    if(check_memory(chi_result))return;
    double estimated_value=info->height*info->width/256;

    for(i=0;i<=255;i++)
    {
        chi_result[0]+=(R[i]-estimated_value)*(R[i]-estimated_value)/estimated_value;
        chi_result[1]+=(G[i]-estimated_value)*(G[i]-estimated_value)/estimated_value;
        chi_result[2]+=(B[i]-estimated_value)*(B[i]-estimated_value)/estimated_value;
    }
    printf("Imaginea %s:\nR:%.2f\nG:%.2f\nB:%.2f\n",type,chi_result[0],chi_result[1],chi_result[2]);
    free(image);free(R);free(G);free(B);free(info);free(chi_result);
    fclose(in);

}
pixel **get_matrix(FILE *in,BMP_data *info)
{
    pixel **matrix=(pixel**)malloc(info->height*sizeof(pixel*));
    if(check_memory(matrix))return NULL;
    int i,j;
    for(i=0;i<info->height;i++)
    {
        matrix[i]=(pixel *)malloc(info->width*sizeof(pixel));
        if(check_memory(matrix[i]))return NULL;
    }

    fseek(in,54,SEEK_SET);
    for (i = info->height - 1; i >= 0; i--)
    {
        for(j=0;j<info->width;j++)
        {
            fread(&matrix[i][j].B,1,1,in);
            fread(&matrix[i][j].G,1,1,in);
            fread(&matrix[i][j].R,1,1,in);
        }
        fseek(in,info->padding,SEEK_CUR);
        if(i==0)
            break;
    }
    return matrix;
}
void load(pixel **image, BMP_data *info, FILE *out)
{

    fwrite(info->header,54,1,out);
    int i, j;
    unsigned char a=0;
    for(i=info->height-1;i>=0;i--)
    {
        for(j=0;j<info->width;j++)
        {
            fwrite(&image[i][j].B,1,1,out);
            fwrite(&image[i][j].G,1,1,out);
            fwrite(&image[i][j].R,1,1,out);
        }
        fwrite(&a,1,info->padding,out);
    }
}
pixel **grayscale(pixel **initial_matrix, BMP_data *info)
{
    pixel **new_matrix=(pixel**)malloc(info->height*sizeof(pixel*));
    if(check_memory(new_matrix))return NULL;
    int i,j;
    for(i=0;i<info->height;i++)
    {
        new_matrix[i]=(pixel *)malloc(info->width*sizeof(pixel));
        if(check_memory(new_matrix[i]))return NULL;
    }
    for(i=0;i<info->height;i++)
    {
        for(j=0;j<info->width;j++)
        {
            new_matrix[i][j].B=new_matrix[i][j].G=new_matrix[i][j].R=
            0.299*initial_matrix[i][j].R+0.587*initial_matrix[i][j].G+
            0.114*initial_matrix[i][j].B;
        }
    }
    return new_matrix;
}
double average_templ(pixel **matrix, BMP_data *info)
{
    double val=0;
    int i, j;
    for(i=0;i<info->height;i++)
    {
        for(j=0;j<info->width;j++)
        {
            val+=matrix[i][j].R;
        }
    }
    val=val/(info->height*info->width);
    return val;
}
double average(pixel **image,int x, int y, BMP_data *info_image, BMP_data *info_template)
{
    double val=0;
    int i, j;

    for(i=x;i<x+info_template->height;i++)
    {
        if(i<info_image->height)
            {for(j=y;j<y+info_template->width;j++)
            {
                if(j<info_image->width)
                    val+=image[i][j].R;
            }}
    }
    val=val/(info_template->height*info_template->width);
    return val;
}
double standard_deviation_templ(pixel **matrix, BMP_data *info, double average)
{
    double deviation, sum=0;
    int i, j;

    for(i=0;i<info->height;i++)
    {
        for(j=0;j<info->width;j++)
        {
            sum+=(matrix[i][j].R-average)*(matrix[i][j].R-average);
        }
    }
    deviation=sqrt(sum/(info->width*info->height-1));

    return deviation;
}
double standard_deviation(pixel **image, int x, int y, BMP_data *info_image, BMP_data *info_template, double average)
{
    double deviation,sum=0;
    int i,j;
    for(i=x;i<x+info_template->height;i++)
    {
        if(i<info_image->height)
            for(j=y;j<y+info_template->width;j++)
        {
            if(j<info_image->width)
                sum+=(image[i][j].R-average)*(image[i][j].R-average);
            else sum+=average*average;
        }
        else sum+=average*average;
    }
    deviation=sqrt(sum/(info_template->width*info_template->height-1));

    return deviation;

}
double correlation(pixel **image,int x,int y,pixel **templ,BMP_data *info_image,BMP_data *info_template,double average_template,double deviation_template)
{
    double sum=0, average_image;
    double deviation_image;
    int i,j;

    average_image=average(image,x,y,info_image,info_template);
    deviation_image=standard_deviation(image,x,y,info_image,info_template,average_image);

    for(i=x;i<x+info_template->height;i++)
    {
        if(i<info_image->height)
        {
            for(j=y;j<y+info_template->width;j++)
            {
                if(j<info_image->width)
                    sum+=(image[i][j].R-average_image)*(templ[i-x][j-y].R-average_template);
                else sum+=(-average_image)*(templ[i-x][j-y].R-average_template);
            }
        }
        else
            for(j=0;j<info_template->width;j++)
            sum+=(-average_image)*(templ[i-x][j].R-average_template);
    }
    double corr;
    corr=sum/(info_template->height*info_template->width*deviation_image*deviation_template);
    return corr;
}
pixel *init_color()
{
    pixel *color=(pixel*)malloc(sizeof(pixel)*10);
    if(check_memory(color))return NULL;
    color[0].R=255; color[0].G=0; color[0].B=0;
    color[1].R=255; color[1].G=255; color[1].B=0;
    color[2].R=0; color[2].G=255; color[2].B=0;
    color[3].R=0; color[3].G=255; color[3].B=255;
    color[4].R=255; color[4].G=0; color[4].B=255;
    color[5].R=0; color[5].G=0; color[5].B=255;
    color[6].R=192; color[6].G=192; color[6].B=192;
    color[7].R=255; color[7].G=140; color[7].B=0;
    color[8].R=128; color[8].G=0; color[8].B=128;
    color[9].R=128; color[9].G=0; color[9].B=0;
    return color;
}
void draw(pixel **image, int x,int y, pixel color, BMP_data *info_image, BMP_data *info_template)
{
    ///se coloreaza prima linie
    for(int j=y;j<y+info_template->width;j++)
        if(x<info_image->height&&j<info_image->width)
            image[x][j]=color;
    ///se coloreaza ultima linie
    for(int j=y;j<y+info_template->width;j++)
        if(x+info_template->height-1<info_image->height&&j<info_image->width)
            image[x+info_template->height-1][j]=color;
    ///se coloreaza prima coloana
    for(int i=x;i<x+info_template->height;i++)
        if(y<info_image->width&&i<info_image->height)
            image[i][y]=color;
    ///se coloreaza ultima coloana
    for(int i=x;i<x+info_template->height;i++)
        if(y+info_template->width-1<info_image->width&&i<info_image->height)
            image[i][y+info_template->width-1]=color;
}
window_list *one_templ_maching(pixel **image, pixel **templ, double ps, BMP_data *info_image, BMP_data *info_template,int digit)
{
    window_list *detections=(window_list*)malloc(info_image->height*info_image->width*sizeof(window_list));
    if(check_memory(detections))return NULL;
    ///Transformarea sablonului -> grayscale
    templ=grayscale(templ,info_template);
    ///Calcularea deviatiei standard si a mediei sablonului
    int n=0;
    double correl,avg,dev;
    avg=average_templ(templ,info_template);
    dev=standard_deviation_templ(templ,info_template,avg);
    ///Centrarea si calcularea corelatiei
    for(int i=0;i<info_image->height;i++)
    {
        for(int j=0;j<info_image->width;j++)
        {
            correl=correlation(image,i,j,templ,info_image,info_template,avg,dev);

            if(correl>ps)
            {
                detections[n].corr=correl;
                detections[n].line=i;
                detections[n].col=j;
                detections[n].elim=0;
                detections[n].digit=digit;
                n++;
            }
        }
    }

    detections[0].size=n;
    return detections;
}
int cmp(const void *a,const void *b)
{
    return ((*((window_list*)a)).corr > (*((window_list*)b)).corr)? -1 : 1;
}
int check_overlap(int x1, int y1, int x2, int y2, BMP_data info1, BMP_data info2)
{
    if(x1>x2+info2.height||x2>x1+info1.height)
        return 0;
    if(y1>y2+info2.width||y2>y1+info2.height)
        return 0;
    return 1;
}
double calculate_overlap(int x1, int y1,int x2, int y2, BMP_data info1, BMP_data info2)
{
    double length, width;
    if(check_overlap(x1,y1,x2,y2,info1,info2)==0)
        return 0;
    if(y1<=y2)///Prima fereastra e in stanga celei de-a doua
        length=y1+info1.width-y2;
    else///Invers
        length=y2+info2.width-y1;
    if(x1<=x2)///Prima fereastra e deasupra celei de-a doua
        width=x1+info1.height-x2;
    else///Invers
        width=x2+info2.height-x1;

    return (double)(length*width)/(info1.height*info1.width+info2.height*info2.width-length*width);

}
BMP_data **get_data_template(char **templates,int n)
{
    BMP_data **info=(BMP_data**)malloc(n*sizeof(BMP_data*));
    if(check_memory(info))return NULL;
    int i;
    for(i=0;i<n;i++)
    {
        info[i]=(BMP_data*)malloc(sizeof(BMP_data));
        if(check_memory(info[i]))return NULL;
    }
    FILE *in;
    for(i=0;i<n;i++)
    {
        in=fopen(templates[i],"r");
        if(in==NULL)
        {
            printf("ERROR:FILE");
            return NULL;
        }
        info[i]=get_data(in);
        fclose(in);
    }
    return info;

}
void elim_non_maxim(window_list *detections,unsigned long n,pixel **image,BMP_data **info)
{
    double overlap;
    ///Sortare in ordine descrescatoare
    qsort(detections,n,sizeof(window_list),cmp);
    ///Eliminare
    for(unsigned long i=0;i<n-1;i++)
    {
        if(detections[i].elim==0)
        {
             for(unsigned long j=i+1;j<n;j++)
            {
                if(detections[j].elim==0)
                {
                    overlap=calculate_overlap(detections[i].line,detections[i].col,detections[j].line,detections[j].col,*(info[detections[i].digit]),*(info[detections[j].digit]));
                    if(overlap>0.2)
                        detections[j].elim=1;
                }
            }
        }
    }
}
void template_matching(FILE *images_name, char *in_out, int n, char **templates, BMP_data **info_templates)
{
    double ps=0.5;
    int i,j;
    ///Retinerea pixelilor din imagine sub forma de matrice
    FILE *bmp_in=fopen(in_out,"rb");
    if(bmp_in==NULL)
    {
        printf("ERROR:FILE");
        return;
    }
    BMP_data *bmp_info=get_data(bmp_in);

    pixel **bmp_image=(pixel**)malloc(bmp_info->height*sizeof(pixel));
    if(check_memory(bmp_image))return;
    for(i=0;i<bmp_info->height;i++)
    {
        bmp_image[i]=(pixel*)malloc(bmp_info->width*sizeof(pixel));
        if(check_memory(bmp_image[i]))return;
    }
    bmp_image=get_matrix(bmp_in,bmp_info);
    ///Transformare - Grayscale -> imaginea principala
    bmp_image=grayscale(bmp_image,bmp_info);
    ///Retinerea sabloanelor sub forma de matrice
    FILE *template_in;
    pixel ***bmp_template=(pixel***)malloc(n*sizeof(pixel**));
    for(i=0;i<n;i++)
    {
        template_in=fopen(templates[i],"rb");
        if(template_in==NULL)
        {
            printf("ERROR:FILE");
            return;
        }
        bmp_template[i]=(pixel**)malloc(info_templates[i]->height*sizeof(pixel*));
        if(check_memory(bmp_template[i]))return;
        for(j=0;j<info_templates[i]->height;j++)
        {
            bmp_template[i][j]=(pixel*)malloc(info_templates[i]->width*sizeof(pixel));
            if(check_memory(bmp_template[i][j]))return;
        }

        bmp_template[i]=get_matrix(template_in,info_templates[i]);
        fclose(template_in);
    }
    ///Recunoasterea cifrelor de mana
    window_list *detections=(window_list*)malloc(n*bmp_info->width*bmp_info->height*sizeof(window_list));
    window_list *aux=(window_list*)malloc(bmp_info->width*bmp_info->height*sizeof(window_list));
    unsigned long k=0;
    for(int i=0;i<n;i++)
    {
        aux=one_templ_maching(bmp_image,bmp_template[i],ps,bmp_info,info_templates[i],i);
        for(j=k;j<k+aux[0].size;j++)
            detections[j]=aux[j-k];
        k=k+aux[0].size;
    }
    ///Eliminarea non-maximelor
    elim_non_maxim(detections,k,bmp_image,info_templates);
    ///Desenarea dreptunghiurilor colorate
    pixel *color;
    color=init_color();
    bmp_image=get_matrix(bmp_in,bmp_info);
    for(i=0;i<k;i++)
        if(detections[i].elim==0)
            draw(bmp_image,detections[i].line,detections[i].col,color[detections[i].digit],bmp_info,info_templates[detections[i].digit]);
    ///Afisare imagine
    fscanf(images_name,"%s",in_out);
    FILE *out=fopen(in_out,"wb");
    if(out==NULL)
    {
        printf("ERROR:FILE");
        return;
    }
    load(bmp_image,bmp_info,out);
    printf("Gata!");
    ///Inchiderea fisierelor si eliberarea memoriei
    fclose(images_name); fclose(bmp_in); fclose(out);
    free(color); free(in_out); free(templates); free(info_templates);
    free(bmp_info); free(bmp_template); free(bmp_image);
    free(detections); free(aux);
}
int main()
{   ///Criptare/decriptare
    ///Citirea denumirii fisierului in care se regasec caile
    char *file_name=(char*)malloc(200*sizeof(char));
    if(check_memory(file_name))return 0;
    printf("Denumirea fisierului in care se gasesc caile celor 2 imagini(cea de criptat si \ncea criptata) si a fisierului care contine cheia secreta:");
    scanf("%s",file_name);
    fflush(stdin);
    FILE *image_encrypt=fopen(file_name,"r");
    if(image_encrypt==NULL)
    {
        printf("ERROR:FILE");
        return 0;
    }
    ///Citirea cailor
    char *image=(char*)malloc(200*sizeof(char));
    char *crypted=(char*)malloc(200*sizeof(char));
    char *decrypted=(char*)malloc(200*sizeof(char));
    char *secret_key=(char*)malloc(200*sizeof(char));
    if(check_memory(image)||check_memory(crypted)||check_memory(decrypted)||check_memory(secret_key))return 0;
    fscanf(image_encrypt,"%s",image);
    fscanf(image_encrypt,"%s",crypted);
    fscanf(image_encrypt,"%s",secret_key);
    ///Criptare
    encrypt(image,crypted,secret_key);
    ///Decriptare
    printf("Denumirea fisierului in care se gasesc caile celor 2 imagini(cea criptata si \ncea decriptata) si a fisierului care contine cheia secreta:");
    scanf("%s",file_name);
    fflush(stdin);
    FILE*image_decrypt=fopen(file_name,"r");
    if(image_decrypt==NULL)
    {
        printf("ERROR:FILE");
        return 0;
    }
    fscanf(image_decrypt,"%s",crypted);
    fscanf(image_decrypt,"%s",decrypted);
    fscanf(image_decrypt,"%s",secret_key);
    decrypt(crypted,decrypted,secret_key);
    ///Testul chi-patrat
    chi_test(image,"initiala");
    chi_test(crypted,"criptata");
    printf("\n");

    ///Template matching
    ///Citirea denumirii fisierului in care se regasec caile
    printf("Denumirea fisierului in care se regasesc caile imaginii si a sabloanelor:\n");
    scanf("%s",file_name);
    fflush(stdin);
    ///Deschiderea fisierului
    FILE *images_name=fopen(file_name,"r");
    if(images_name==NULL)
    {
        printf("ERROR:FILE");
        return 0;
    }
    char *in_out=(char*)malloc(200*sizeof(char));
    if(check_memory(in_out))return 0;
    fscanf(images_name,"%s",in_out);

    int n,i,j;
    ///Citirea numarului de sabloane
    printf("\nCate sabloane? n=");
    scanf("%d",&n);
    printf("\nDureaza ceva...");
    ///Citirea denumirilor sabloanelor
    char **templates=(char**)malloc(n*sizeof(char*));
    if(check_memory(templates))return 0;
    for(i=0;i<n;i++)
    {
        templates[i]=(char*)malloc(200*sizeof(char));
        if(check_memory(templates[i]))return 0;
        fscanf(images_name,"%s",templates[i]);
    }
    ///Citirea informatiilor despre fiecare sablon
    BMP_data **info_templates;
    info_templates=get_data_template(templates,n);
    ///Rularea functiei de template matching
    template_matching(images_name,in_out,n,templates,info_templates);

    return 0;
}
