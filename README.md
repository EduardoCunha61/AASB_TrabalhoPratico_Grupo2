# AASB_TrabalhoPratico_Grupo2
Relatório do trabalho AASB - Protein Sequence Manager
Ana Barbosa, 
Cátia Silva, 
Eduardo Cunha.

Inicialmente, adaptou-se a classe MySeq das aulas para este trabalho, criando a classe Seq().
Cada objeto tem um atributo sequência(string) e um dicionário com as propriedades referentes a cada ao respetivo objeto.
As propriedades guardadas são as seguintes: "Id NCBI","Id Uniprot","Function","Genes","Taxonomic lineage" e Size da sequencia (calculado automaticamente).
Toda a informação é guardada num dicionário, sendo que a chave para cada objecto Seq é um número inteiro.

A aplicação importa e exporta a informação na "base de dados" para ficheiros .txt, que são definidos no ficheiro main.py como variáveis globais.

Para correr a aplicação, basta efetuar o comando python3 main.py. Também é necessário instalar o módulo Biopython através do
comando pip install Biopython (que será usado nas funcionalidades mais avançadas do projeto).

Este trabalho está definido em três ficheiros python diferentes: 
-> Main.py: Main da aplicação. Todas as operações IO ocorrem neste ficheiro, assim como as inserções e acessos à "base de dados".
-> Functionalities.py: Ficheiro com funções avançadas como  Importar e exportar ficheiros, preparar os objetos para criar um novo objeto seq, etc.
-> Classes.py: Definição de todas as classes que são usadas ao longo do trabalho.

Relativamente aos requisitos obrigatórios, todos estão implementados corretamente.
Já os requisitos extra, para os pontos em que era necessário recorrer a programas externos, como o caso do Alinhamento múltiplo, foi necessária a instalação de um software abordado no decorrer do curso - ClustalW e, para o último ponto, que requeria a procura de modificações pós-tradução, decidimos adaptar a procura não por sequência mas sim por identificador, neste caso em concreto, o identificador da Base de Dados UniProt.
