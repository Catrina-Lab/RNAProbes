FROM alpine/curl

CMD curl https://${REDIS_LOCATION}.upstash.io/incr/auto_updates -H "Authorization: Bearer ${REDIS_AUTH_TOKEN}"

#sudo docker build  --tag update:redis --no-cache .
#sudo /home/linuxbrew/.linuxbrew/bin/flyctl deploy --local-only --build-only --push --dockerfile Dockerfile -a pingredis
#sudo /home/linuxbrew/.linuxbrew/bin/flyctl m run . --name pingredis --config fly.toml --schedule "daily"

#if it doesn't work, try updating secrets in fly (click secrets, deploy to all machines)