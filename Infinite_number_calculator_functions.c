#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

struct Node{
    int* mat;
    int sign;
    int length;
};typedef struct Node node;

node* GetInt(char s[]) {
    int n = strlen(s);
    node* arr= (node*)malloc(sizeof(node));
    char oper[]="-";
    if(s[0] == '-'){
        arr->sign=1;
        int i=0;
        n=strlen(s)-1;
        arr->length=n;
        arr->mat=(int*)malloc(n*sizeof(int));
        while(i<n){
            arr->mat[i]=s[i+1]-'0';
            i++;
        }
    }else {
        arr->mat=(int*)malloc(n*sizeof(int));
        for (int i = 0; i < n; i++) {
            arr->mat[i] = s[i] - '0';
        }
        arr->sign=0;
        arr->length=n;
        //printf("%d",arr->sign);
    }
    return arr;
}

char* GetChar(node* arr){
    if(arr->sign==1){
        char* s = (char *) malloc(arr->length*sizeof(char)+2);
        char oper[]="-";
        s[0]=oper[0];
        int index=1;
        for (int i = 0; i < arr->length; i++) {
            index+=sprintf(&s[index],"%d",arr->mat[i]) ;
        }
        return s;
    }else {
        char* s = (char *) malloc(arr->length * sizeof(char)+1);
        int index = 0;
        for (int i = 0; i < arr->length; i++) {
            index += sprintf(&s[index], "%d", arr->mat[i]);
        }
        return s;
    }
}

char* Add(node* num1,node* num2){
    int m=MAX(num1->length,num2->length);
    int l=MIN(num1->length,num2->length);
    int maloc_size=m+1;
    int k=m-l;
    int add =1;
    if (num1->mat[0]==0 || num2->mat[0]==0){
        maloc_size=m;
        add=0;
    }
    node* res=(node*)calloc(1,sizeof(node));
    res->mat=(int*)malloc(maloc_size*sizeof(int));
    //int res[m+1];
    int i=l-1;
    int helper=0;
    while(i>-1 ){
        if(num2->length==l) {
            res->mat[k + i + add] = num1->mat[k + i] + num2->mat[i] + helper;
        }
        if(num1->length==l){
            res->mat[k + i + add] = num2->mat[k + i] + num1->mat[i] + helper;
        }
        if(res->mat[k+i+add]>9){
            res->mat[k+i+add]=res->mat[k+i+add]%10 ;
            helper=1;
        }else{
            helper=0;
        }
        i--;
    }

    while(k>0) {
        if (num1->length == l) {
            if (add==0)
            {
                res->mat[k-1] = num2->mat[k-1]+helper;
            } else{
                res->mat[k] = num2->mat[k-1]+helper;
            }

        }
        if (num2->length == l) {
            if (add==0){
                res->mat[k-1] = num1->mat[k-1]+helper;
            }else{
                res->mat[k] = num1->mat[k-1]+helper;
            }

        }
        if (res->mat[k] > 9) {
            res->mat[k] =res->mat[k]% 10;
            helper = 1;
        }else{
            helper=0;
        }
        k--;
    }
    res->mat[0]=helper;
  if(res->mat[0]==0) {
       res->mat = res->mat+1;
       maloc_size=m;
   }
  res->length=maloc_size;
   if(num1->sign==1){res->sign=1;}
   if(num1->sign==0){res->sign=0;}
        char *s = GetChar(res);
        free(res->mat);
        free(res);
        return s;

}

char* mul(node* num1,node* num2){
    node* a=(node*)malloc(sizeof(node));
    a->mat=(int*)calloc(num1->length,sizeof(int)) ;
   // if(!a){printf("Error mem alloc for a"); }
   node* b=(node*)malloc(num2->length*sizeof(node));
   b->mat =(int*)calloc(num2->length,sizeof(int)) ;
   // if(!b){printf("Error mem alloc for b"); }
   node* r=(node*)calloc(1,sizeof(node));
   r->mat=(int*)calloc((num1->length+num2->length),sizeof(int)) ;
   // if(!r){printf("Error mem alloc for r"); }
    for(int k=num1->length-1,m=0;k>-1;k--,m++){
        a->mat[m]=num1->mat[k];
    }
    for(int k=num2->length-1,m=0;k>-1;k--,m++){
        b->mat[m]=num2->mat[k];
    }
    node* ans=(node*)calloc(1,sizeof(node));
    ans->mat=(int*)calloc((num1->length+num2->length),sizeof(int));
    //if(!ans){printf("Error mem alloc for ans"); }
    for(int i=0;i<num1->length;i++){
        for(int j=0;j<num2->length;j++){
            ans->mat[i+j]+=a->mat[i]*b->mat[j];
        }
    }
    for(int i=0;i<num1->length+num2->length-1;i++){
        int temp=ans->mat[i]/10;
        ans->mat[i]=ans->mat[i]%10;
        ans->mat[i+1]=ans->mat[i+1]+ temp;
    }
    for(int i=num1->length+num2->length-1,l=0;i>-1;i--,l++){
        r->mat[l]=ans->mat[i];
    }
    free(a->mat);
    free(a);
    free(b->mat);
    free(b);
    free(ans->mat);
    free(ans);
    r->length=num1->length+num2->length;
    if(num1->sign!=num2->sign){r->sign=1;}
    if(num1->sign==num2->sign){r->sign=0;}

  if(r->mat[0]==0){
      r->mat=r->mat+1;
    r->length=num1->length+num2->length-1;
   char* st=GetChar(r);
   free(r->mat);
   free(r);
   return st;
  }else {
      char *str = GetChar(r);
      free(r->mat);
      free(r);
      return str;
  }
}

int g(node* num1,node* num2){
    if(num1->sign>num2->sign){return 0;}
    if(num1->sign<num2->sign){return 1;}
    if(num1->length>num2->length){return 1;}
    if(num1->length<num2->length){return 0;}
    else{
        for(int i=0;i<num1->length;i++){
            if(num1->mat[i]>num2->mat[i]){return 1;}
            if(num2->mat[i]>num1->mat[i]){return 0;}
        }
        printf("There are equal");
    }
}

int gabs(node* num1,node* num2){
    if(num1->length>num2->length){return 1;}
    if(num1->length<num2->length){return 0;}
    else{
        for(int i=0;i<num1->length;i++){
            if(num1->mat[i]>num2->mat[i]){return 1;}
            if(num2->mat[i]>num1->mat[i]){return 0;}
        }
        printf("There are equal");
    }
}

int e(node* num1,node* num2){
    if(num1->length!=num2->length){return 0;}
    if(num1->sign!=num2->sign){return 0;}
    int i=0;
    while(num1->mat[i]==num2->mat[i] && i<num1->length){
        i++;
    }
    if(i!=num1->length){return 0;}
    else {return 1;}
}

char* sub(node* num1,node* num2){
    int m=MAX(num1->length,num2->length);
    int l=MIN(num1->length,num2->length);
    int k=m-l;
    node* res=(node*)malloc(sizeof(node));
    res->mat=(int*)malloc(m*sizeof(int));
    int helper=0;
    if(m==num1->length && (num1->length!=num2->length)){
        int i=l-1;
        while(i>-1) {
            if (num1->mat[k + i] < (num2->mat[i] + helper)) {
                res->mat[k + i] = num1->mat[k + i] + 10 - (num2->mat[i] + helper);
                helper = 1;
            } else {
                res->mat[k + i] = num1->mat[k + i] - (num2->mat[i] + helper);
                helper = 0;
            }
            i--;
        }
        while(k>0) {
            res->mat[k - 1] = num1->mat[k - 1]-helper;
            helper=0;
            k--;
        }
        res->sign=num1->sign;
        res->length=m;
        if(res->mat[0]==0){
            res->mat=res->mat+1;
            res->length=m-1;
        }
            char *s = GetChar(res);
        free(res->mat);
        free(res);
            return s;

    }
////////////////////////////////////////////////////////
    if(m==num2->length && num1->length!=num2->length) {
        int i = l - 1;
        while (i > -1) {
            if (num2->mat[k + i] < (num1->mat[i] + helper)) {
                res->mat[k + i] = num2->mat[k + i] + 10 - (num1->mat[i] + helper);
                helper = 1;
            } else {
                res->mat[k + i] = num2->mat[k + i] - (num1->mat[i] + helper);
                helper = 0;
            }
            i--;
        }
        while (k > 0) {
            res->mat[k - 1] = num2->mat[k - 1] - helper;
            helper = 0;
            k--;
        }

        res->sign=num2->sign;
        res->length=m;
        if (res->mat[0] == 0) {
            res->mat=res->mat+1;
            res->length=m-1;
            char *s = GetChar(res);
            free(res->mat);
            free(res);
            return s;
        }
            char *s = GetChar(res);
        free(res->mat);
        free(res);
            return s;

    }
////////////////////////////////////////////////////
    if(num1->length==num2->length){
        if(e(num1,num2)==1){
            char*s =NULL;//(char*)malloc(0*sizeof(int));
            return s;
        }
        if(gabs(num1,num2)==1){
            res->sign=num1->sign;
            int i=num1->length-1;
            helper=0;
            while(i>-1) {
                if (num1->mat[i] < (num2->mat[i] + helper)) {
                    res->mat[i] = num1->mat[i] + 10 - (num2->mat[i] + helper);
                    helper = 1;
                } else {
                    res->mat[i] = num1->mat[i] - (num2->mat[i] + helper);
                    helper = 0;
                }
                i--;
            }
            int c=0; int b=0;
        while(res->mat[b]==0){
            c=c+1;
            b++;
        }
        if(c) {
            printf("%d\n",c);
            res->mat=res->mat+c;
            res->length=m-c;
            char* s=GetChar(res);
            free(res->mat);
            free(res);
            return s;
        } else{
            res->length=m;
            char* s=GetChar(res);
            free(res->mat);
            free(res);
            return s;
        }
        }else{
            res->sign=num2->sign;
            res->length=m;
            int i=num1->length-1;
            helper=0;
            while (i > -1) {
                if (num2->mat[i] < (num1->mat[i] + helper)) {
                    res->mat[i] = num2->mat[i] + 10 - (num1->mat[i] + helper);
                    helper = 1;
                } else {
                    res->mat[i] = num2->mat[i] - (num1->mat[i] + helper);
                    helper = 0;
                }
                i--;
            }
            int c=0;
            while(res->mat[c]==0){
                c=c+1;
            }
            if(c) {
                res->mat = res->mat + c;
                res->length = m - c;

                char *s = GetChar(res);
                free(res->mat);
                free(res);
                return s;
            }
             else{
                char* s=GetChar(res);
                free(res->mat);
                free(res);
                return s;
            }
        }

    }

}

int main() {
    char* s="10";
    node* res=GetInt(s);
    char* s2="10";
    node* res2=GetInt(s2);
// char* str=Add(res,res2);
 //  char* str = mul(res,res2);
 //  char* str=sub(res,res2);
//  int str= e(res,res2);
//  int str=g(res,res2);
char *str=0;
char oper='m';
int t=0;
switch(oper) {
    case 'a':
        if (res->sign == res2->sign) {
            str = Add(res, res2);
            break;
        }
        else{
            str = sub(res, res2);
            break;
        }
    case 'm':
        str=mul(res,res2);
        break;
    case 'e':
        t=e(res,res2);
        break;

    case 'g':
        t=g(res,res2);
        break;

    case 's':
        if(res2->sign != res->sign) {
            str = Add(res, res2);
            break;
        }else {
            if(res2->sign==1)
            {res2->sign=0;}
            else {res2->sign=1;}
            str=sub(res,res2);
            break;
        }
}

printf("%s",str);
//printf("%d",t);
    return 0;
}


