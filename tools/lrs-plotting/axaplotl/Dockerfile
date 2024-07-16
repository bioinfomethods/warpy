FROM node:latest

RUN mkdir /app
WORKDIR /app

ENV SSL_CERT_FILE="/usr/local/share/ca-certificates/ZscalerRootCertificate-2048-SHA256.crt"
RUN curl -o "$SSL_CERT_FILE" http://repo.mcri.edu.au/ZCC/Certs/ZscalerRootCerts/ZscalerRootCertificate-2048-SHA256.crt
ENV NODE_EXTRA_CA_CERTS="$SSL_CERT_FILE"

RUN npm install -g vite
RUN npm install -g serve

COPY . .

RUN npm install
RUN npm run build

CMD ["npm", "run", "preview"]
