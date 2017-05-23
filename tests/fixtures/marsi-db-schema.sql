--
-- PostgreSQL database dump
--

SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SET check_function_bodies = false;
SET client_min_messages = warning;

--
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


SET search_path = public, pg_catalog;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: alembic_version; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE alembic_version (
    version_num character varying(32) NOT NULL
);


ALTER TABLE public.alembic_version OWNER TO postgres;

--
-- Name: metabolite_fingerprints; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE metabolite_fingerprints (
    id integer NOT NULL,
    metabolite_id integer,
    fingerprint_type character varying(10) NOT NULL,
    fingerprint character varying(2048) NOT NULL
);


ALTER TABLE public.metabolite_fingerprints OWNER TO postgres;

--
-- Name: metabolite_fingerprints_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE metabolite_fingerprints_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.metabolite_fingerprints_id_seq OWNER TO postgres;

--
-- Name: metabolite_fingerprints_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE metabolite_fingerprints_id_seq OWNED BY metabolite_fingerprints.id;


--
-- Name: metabolite_references; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE metabolite_references (
    metabolite_id integer,
    reference_id integer
);


ALTER TABLE public.metabolite_references OWNER TO postgres;

--
-- Name: metabolite_synonyms; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE metabolite_synonyms (
    metabolite_id integer,
    synonym_id integer
);


ALTER TABLE public.metabolite_synonyms OWNER TO postgres;

--
-- Name: metabolites; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE metabolites (
    id integer NOT NULL,
    inchi_key character varying(27) NOT NULL,
    inchi character varying(5000) NOT NULL,
    analog boolean,
    formula character varying(500) NOT NULL,
    num_atoms integer NOT NULL,
    num_bonds integer NOT NULL,
    sdf text,
    solubility double precision,
    num_rings integer NOT NULL
);


ALTER TABLE public.metabolites OWNER TO postgres;

--
-- Name: metabolites_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE metabolites_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.metabolites_id_seq OWNER TO postgres;

--
-- Name: metabolites_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE metabolites_id_seq OWNED BY metabolites.id;


--
-- Name: references; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE "references" (
    id integer NOT NULL,
    database character varying(100) NOT NULL,
    accession character varying(100) NOT NULL
);


ALTER TABLE public."references" OWNER TO postgres;

--
-- Name: references_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE references_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.references_id_seq OWNER TO postgres;

--
-- Name: references_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE references_id_seq OWNED BY "references".id;


--
-- Name: synonyms; Type: TABLE; Schema: public; Owner: postgres; Tablespace: 
--

CREATE TABLE synonyms (
    id integer NOT NULL,
    synonym character varying(500) NOT NULL
);


ALTER TABLE public.synonyms OWNER TO postgres;

--
-- Name: synonyms_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE synonyms_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.synonyms_id_seq OWNER TO postgres;

--
-- Name: synonyms_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE synonyms_id_seq OWNED BY synonyms.id;


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY metabolite_fingerprints ALTER COLUMN id SET DEFAULT nextval('metabolite_fingerprints_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY metabolites ALTER COLUMN id SET DEFAULT nextval('metabolites_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY "references" ALTER COLUMN id SET DEFAULT nextval('references_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY synonyms ALTER COLUMN id SET DEFAULT nextval('synonyms_id_seq'::regclass);


--
-- Name: _database_accession_uc; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY "references"
    ADD CONSTRAINT _database_accession_uc UNIQUE (database, accession);


--
-- Name: _fp_type_uc; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY metabolite_fingerprints
    ADD CONSTRAINT _fp_type_uc UNIQUE (metabolite_id, fingerprint_type);


--
-- Name: alembic_version_pkc; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY alembic_version
    ADD CONSTRAINT alembic_version_pkc PRIMARY KEY (version_num);


--
-- Name: metabolite_fingerprints_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY metabolite_fingerprints
    ADD CONSTRAINT metabolite_fingerprints_pkey PRIMARY KEY (id);


--
-- Name: metabolites_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY metabolites
    ADD CONSTRAINT metabolites_pkey PRIMARY KEY (id);


--
-- Name: references_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY "references"
    ADD CONSTRAINT references_pkey PRIMARY KEY (id);


--
-- Name: synonyms_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres; Tablespace: 
--

ALTER TABLE ONLY synonyms
    ADD CONSTRAINT synonyms_pkey PRIMARY KEY (id);


--
-- Name: uq_inchi_key; Type: INDEX; Schema: public; Owner: postgres; Tablespace: 
--

CREATE UNIQUE INDEX uq_inchi_key ON metabolites USING btree (inchi_key);


--
-- Name: uq_synonyms; Type: INDEX; Schema: public; Owner: postgres; Tablespace: 
--

CREATE UNIQUE INDEX uq_synonyms ON synonyms USING btree (synonym);


--
-- Name: metabolite_fingerprints_metabolite_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY metabolite_fingerprints
    ADD CONSTRAINT metabolite_fingerprints_metabolite_id_fkey FOREIGN KEY (metabolite_id) REFERENCES metabolites(id);


--
-- Name: metabolite_references_metabolite_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY metabolite_references
    ADD CONSTRAINT metabolite_references_metabolite_id_fkey FOREIGN KEY (metabolite_id) REFERENCES metabolites(id);


--
-- Name: metabolite_references_reference_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY metabolite_references
    ADD CONSTRAINT metabolite_references_reference_id_fkey FOREIGN KEY (reference_id) REFERENCES "references"(id);


--
-- Name: metabolite_synonyms_metabolite_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY metabolite_synonyms
    ADD CONSTRAINT metabolite_synonyms_metabolite_id_fkey FOREIGN KEY (metabolite_id) REFERENCES metabolites(id);


--
-- Name: metabolite_synonyms_synonym_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY metabolite_synonyms
    ADD CONSTRAINT metabolite_synonyms_synonym_id_fkey FOREIGN KEY (synonym_id) REFERENCES synonyms(id);


--
-- Name: public; Type: ACL; Schema: -; Owner: postgres
--

REVOKE ALL ON SCHEMA public FROM PUBLIC;
REVOKE ALL ON SCHEMA public FROM postgres;
GRANT ALL ON SCHEMA public TO postgres;
GRANT ALL ON SCHEMA public TO PUBLIC;


--
-- PostgreSQL database dump complete
--

